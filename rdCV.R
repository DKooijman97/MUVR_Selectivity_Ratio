#' rdCV: Wrapper for repeated double cross-validation without variable selection
#'
#' @param X Independent variables. NB: Variables (columns) must have names/unique identifiers. NAs not allowed in data. For ML, X is upper half only (X1-X2)
#' @param Y Response vector (Dependent variable). For DA (classification), Y should be factor or character. For ML, Y is omitted. For regression, Y is numeric.
#' @param ID Subject identifier (for sampling by subject; Assumption of independence if not specified)
#' @param nRep Number of repetitions of double CV.
#' @param nOuter Number of outer CV loop segments.
#' @param nInner Number of inner CV loop segments.
#' @param DA Logical for Classification (discriminant analysis) (Defaults do FALSE, i.e. regression). PLS is limited to two-class problems (see `Y` above).
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'MISS' for classification; or 'RMSEP' (default) for regression.)
#' @param method Multivariate method. Supports 'PLS' and 'RF' (default)
#' @param methParam List with parameter settings for specified MV method (defaults to ???)
#' @param ML Logical for multilevel analysis (defaults to FALSE)
#' @param modReturn Logical for returning outer segment models (defaults to FALSE)
#' @param logg Logical for whether to sink model progressions to `log.txt`
#'
#' @return An object containing stuff...
#' @export
rdCV=function(X,Y,ID,nRep=5,nOuter=6,nInner,DA=FALSE,fitness=c('AUROC','MISS','RMSEP'),method=c('PLS','RF'),methParam,ML=FALSE,modReturn=FALSE,logg=FALSE){
  ### Code adapted from MVWrap - Brokenness and oddities are most likely due to this
  library(pROC)
  # Initialise modelReturn with function call
  modelReturn=list(call=match.call())
  # Start timer
  start.time=proc.time()[3]
  # Check indata
  if (is.null(dim(X))){
    cat('\nError: Wrong format of X matrix.\n')
    return(NULL)
  }
  nSamp=nrow(X)
  nVar=nVar0=ncol(X)
  if (missing(ID)) {
    cat('\nMissing ID -> Assume all unique (i.e. sample independence)')
    ID=1:nSamp
  }
  if (missing(nInner)) nInner=nOuter-1
  if (missing(method)) method='RF'
  if (method=='RF') library(randomForest) else library(mixOmics)
  if (missing(methParam)) {
    if (method=='PLS') {
      methParam=list(compMax=ifelse(nVar<5,nVar,5),mode='regression')
    } else {
      methParam=list(ntreeIn=150,ntreeOut=300,mtryMaxIn=150)
    }
  }
  if (ML) {
    X=rbind(X,-X)
    Y=rep(c(-1,1),each=nSamp)
    nSamp=2*nSamp
    ID=c(ID,ID)
    DA=FALSE
    fitness='MISS'
    cat('\nMultilevel -> Regression on (-1,1) & fitness=MISS')
  }
  if (!is.null(dim(Y))) {
    cat('\nY is not a vector: Return NULL')
    return(NULL)
  }
  if (is.character(Y)) Y=factor(Y)
  if (is.factor(Y)) {
    cat('\nY is factor -> Classification (',length(unique(Y)),' classes)',sep='')
    DA=TRUE
  }
  if (is.numeric(Y) & DA) {
    Y=as.factor(Y)
    cat('\nDA=TRUE -> Y as factor -> Classification (',length(unique(Y)),' classes)',sep='')
  }
  if (missing(fitness)) {
    if (DA) {
      if (length(unique(Y))>2) {
        fitness='MISS'
        cat('\nMissing fitness -> MISS')
      } else {
        fitness='AUROC'
        cat('\nMissing fitness -> AUROC')
      }
    } else {
      fitness='RMSEP'
      cat('\nMissing fitness -> RMSEP')
    }
  }
  if (nrow(X)!=length(Y)) {
    cat('\nMust have same nSamp in X and Y: Return NULL')
    return(NULL)
  }
  ## Store indata in list for later model return
  InData=list(X=X,Y=Y,ID=ID,nRep=nRep,nOuter=nOuter,nInner=nInner,DA=DA,fitness=fitness,method=method,methParam=methParam,ML=ML)
  ## Sort sampling based in subjects and not index
  unik=!duplicated(ID)  # boolean of unique IDs
  unikID=ID[unik]  
  if (DA) {
    unikY=Y[unik]  # Counterintuitive, but needed for groupings by Ynames
    Ynames=sort(unique(Y))  # Find groups
    groups=length(Ynames) # Number of groups
    groupID=list()  # Allocate list for indices of groups
    for (g in 1:groups) { 
      groupID[[g]]=unikID[unikY==Ynames[g]]  # Find indices per group
    }
    yPred=array(dim=c(length(Y),length(levels(Y)),nRep))
    colnames(yPred)=levels(Y)
    dimnames(yPred)[[3]]=paste('Rep',1:nRep,sep='')
    yPredR=matrix(nrow=length(Y),ncol=length(levels(Y)))
    colnames(yPredR)=levels(Y)
  } else {
    yPred=matrix(nrow=length(Y),ncol=nRep)
    colnames(yPred)=paste('Rep',1:nRep,sep='')
    yPredR=numeric(length(Y))
  }
  # Allocate response vectors and matrices for var's, nComp and VIP ranks over repetitions
  nCompRep=missRep=numeric(nRep)
  names(nCompRep)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
  VIPRep=matrix(data=nVar0,nrow=nVar0,ncol=nRep)
  rownames(VIPRep)=colnames(X)
  colnames(VIPRep)=paste(rep('rep',nRep),1:nRep,sep='')
  VAL=matrix(nrow=nOuter,ncol=nRep)
  rownames(VAL)=paste('outSeg',1:nOuter,paste='')
  colnames(VAL)=paste('rep',1:nRep,paste='')
  ## Choose package/core algorithm according to chosen method
  packs=c(ifelse(method=='PLS','mixOmics','randomForest'),'pROC')
  exports=c(ifelse(method=='PLS','plsInner','rfInner'),'vectSamp')
  ## Start repetitions
  # reps=list()
  # for (r in 1:nRep){
  reps=foreach(r=1:nRep, .packages=packs, .export=exports) %dopar% {
    # r=1
    # r=r+1
    if (logg) sink('log.txt',append=TRUE)
    if (modReturn) outMod=list()
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
    if (DA) {
      groupTest=list()  ## Allocate list for samples within group
      for (gT in 1:groups) { 
        groupTest[[gT]]=vectSamp(groupID[[gT]],n=nOuter)  # Draw random samples within group
      }
      allTest=groupTest[[1]] # Add 1st groups to 'Master' sample of all groups
      for (gT in 2:groups) {  # Add subsequent groups
        allTest=allTest[order(sapply(allTest,length))]
        for (aT in 1:nOuter) {
          allTest[[aT]]=sort(c(allTest[[aT]],groupTest[[gT]][[aT]]))
        }
      }
    } else {
      allTest=vectSamp(unikID,n=nOuter)
    }
    nCompOut=numeric(nOuter)
    names(nCompOut)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    VIPOut=matrix(data=nVar0,nrow=nVar0,ncol=nOuter)
    rownames(VIPOut)=colnames(X)
    colnames(VIPOut)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    VALRep=matrix(nrow=nOuter,ncol=1)
    ## Perform outer loop segments -> one "majority vote" MV model per segment
    for (i in 1:nOuter) {   
      # i=1
      # i=i+1
      cat('\n Segment ',i,':',sep='') # Counter
      ## Draw out test set
      testID=allTest[[i]] # Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex=ID%in%testID # Logical for samples included
      xTest=X[testIndex,]
      yTest=Y[testIndex]
      inID=unikID[!unikID%in%testID]  # IDs not in test set
      if (DA) inY=unikY[!unikID%in%testID]  # Counterintuitive, but needed for grouping by Ynames
      ## Allocate variables for later use
      missIn=aucIn=rmsepIn=PRESSIn=nCompIn=matrix(nrow=nInner,ncol=1)
      rownames(rmsepIn)=rownames(PRESSIn)=rownames(missIn)=rownames(aucIn)=rownames(nCompIn)=paste(rep('inSeg',nInner),1:nInner,sep='')
      colnames(rmsepIn)=colnames(PRESSIn)=colnames(missIn)=colnames(aucIn)=colnames(nCompIn)=nVar
      VIPInner=matrix(data=nVar0,nrow=nVar0,ncol=nInner)
      rownames(VIPInner)=colnames(X)
      colnames(VIPInner)=paste(rep('inSeg',nInner),1:nInner,sep='')
      ## Perform steps with successively fewer variables
        if (method=='PLS') comp=ifelse(nVar<methParam$compMax,nVar,methParam$compMax)
        if (method=='RF') {
          mtryIn=ifelse(DA,
                        ifelse(sqrt(nVar)>methParam$mtryMaxIn,methParam$mtryMaxIn,floor(sqrt(nVar))),
                        ifelse((nVar/3)>methParam$mtryMaxIn,methParam$mtryMaxIn,floor(nVar/3)))
          mtryIn=ifelse(mtryIn<2,2,mtryIn)
        }
        if (DA) {
          groupIDVal=list()
          for (g in 1:groups) { 
            groupIDVal[[g]]=inID[inY==Ynames[g]]  # Find indices per group
          }
          groupVal=list()  ## Allocate list for samples within group
          for (gV in 1:groups) { 
            groupVal[[gV]]=vectSamp(groupIDVal[[gV]],n=nInner)  # Draw random samples within group
          }
          allVal=groupVal[[1]] # Add 1st groups to 'Master' sample of all groups
          for (gV in 2:groups) {  # Add subsequent groups
            allVal=allVal[order(sapply(allVal,length))]
            for (aV in 1:nInner) {
              allVal[[aV]]=sort(c(allVal[[aV]],groupVal[[gV]][[aV]]))
            }
          }
        } else {
          allVal=vectSamp(inID,n=nInner)
        }
        ## Inner CV loop
        for (j in 1:nInner) {
          # j=1 
          # j=j+1
          cat('.') # Counter
          valID=allVal[[j]] # Draw out segment = validation set
          valIndex=ID%in%valID
          xVal=X[valIndex,]
          yVal=Y[valIndex]
          trainID=inID[!inID%in%valID]
          trainIndex=ID%in%trainID # Define Training segment
          xTrain=X[trainIndex,]
          yTrain=Y[trainIndex]
          # sum(trainIndex,valIndex,testIndex)
          # trainIndex|valIndex|testIndex
          ## Make inner model
          if (method=='PLS') {
            inMod=plsInner(xTrain,yTrain,xVal,yVal,DA,fitness,comp,methParam$mode)
            nCompIn[j,1]=inMod$nComp
          } else {
            inMod=rfInner(xTrain,yTrain,xVal,yVal,DA,fitness,ntree=methParam$ntreeIn,mtry=mtryIn)
          }
          # Store fitness metric
          if (fitness=='MISS') {
            missIn[j,1]=inMod$miss
          } else if (fitness=='AUROC') {
            aucIn[j,1]=inMod$auc
          } else {
            rmsepIn[j,1]=inMod$rmsep
            PRESSIn[j,1]=(inMod$rmsep^2)*length(yVal)
          }
          # Store VIPs
          VIPInner[match(names(inMod$vip),rownames(VIPInner)),j]=inMod$vip
      }
      if (fitness=='AUROC') {
        fitRank=-aucIn
        fitRank[]=rank(fitRank)
        fitRank=colMeans(fitRank)
        VALRep[i,]=colMeans(aucIn)
      } else if (fitness=='MISS') {
        fitRank=missIn
        fitRank[]=rank(fitRank)
        fitRank=colMeans(fitRank)
        VALRep[i,]=colSums(missIn)
      }else {
        fitRank=rmsepIn
        fitRank[]=rank(fitRank)
        fitRank=colMeans(fitRank)
        VALRep[i,]=sqrt(colSums(PRESSIn)/sum(!testIndex))
      }
      # Per outer segment: Average inner loop variables, nComp and VIP ranks 
      if (method=='PLS') {
        nCompOut[i]=round(mean(nCompIn[,1]))
      }
      VIPOut[,i]=apply(VIPInner,1,mean)
      # Build outer model for min and max nComp and predict YTEST
      xIn=X[!testIndex,] # Perform Validation on all samples except holdout set
      yIn=Y[!testIndex]
      if (method=='PLS'){
        if (DA) plsOut=plsda(xIn,yIn,ncomp=nCompOut[i]) else 
          plsOut=pls(xIn,yIn,ncomp=nCompOut[i],mode=methParam$mode)
        if (length(plsOut$nzv$Position)>0) removeVar=rownames(plsOut$nzv$Metrics) else removeVar=NA
        incVar=colnames(X)[!colnames(X)%in%removeVar]
        xTest=subset(xTest,select=incVar)
        yPredR[testIndex]=predict(plsOut,newdata=xTest)$predict[,,nCompOut[i]]  # 
        # Prediction of newdata
        if (modReturn) outMod[[i]]=plsOut
      } else {
        rfOut=randomForest(xIn,yIn,xTest,yTest)
        if (DA) {
          yPredR[testIndex,]=rfOut$test$votes
        } else {
          yPredR[testIndex]=rfOut$test$predicted
        }
        if (modReturn) outMod[[i]]=rfOut
      }
    }
    # Per repetition: Average outer loop variables, nComp and VIP ranks 
    parReturn=list(yPred=yPredR)
    parReturn$VIPRep=apply(VIPOut,1,mean)
    if (method=='PLS'){
      parReturn$nCompRep=round(mean(nCompOut))
    }
    parReturn$VAL=VALRep
    if (modReturn) parReturn$outModel=outMod
    if (logg) sink()
    return(parReturn)
    # reps[[r]]=parReturn
  }
  if (modReturn) outMods=list()
  for (r in 1:nRep) {
    if (DA) yPred[,,r]=reps[[r]]$yPred else
      yPred[,r]=reps[[r]]$yPred
    VIPRep[,r]=reps[[r]]$VIPRep
    if (method=='PLS') nCompRep[r]=reps[[r]]$nCompRep
    VAL[,r]=reps[[r]]$VAL
    if (modReturn) outMods=c(outMods,reps[[r]]$outModel)
  }
  # Average predictions
  if (DA) {
    yPredAve=apply(yPred,c(1,2),mean)
  } else {
    yPredAve=apply(yPred,1,mean)
  }
  modelReturn$yPred=yPredAve
  if (DA) {
    auc=numeric(length(levels(Y)))
    names(auc)=levels(Y)
    for (cl in 1:length(levels(Y))) {
      auc[cl]=roc(Y==(levels(Y)[cl]),yPredAve[,cl])$auc
    }
    # Classify predictions
    yClass=as.factor(apply(yPredAve,1,function(x) levels(Y)[which.max(x)]))
    miss=sum(yClass!=Y)
    modelReturn$yClass=yClass
    modelReturn$miss=miss
    modelReturn$auc=auc
  } else if (ML) {
    modelReturn$yClass=ifelse(yPredAve>=0,1,-1)
    modelReturn$miss=sum(modelReturn$yClass!=Y)
    modelReturn$auc=roc(Y,yPredAve)$auc
  }
  # Average VIP ranks over repetitions
  VIP=apply(VIPRep,1,mean)
  modelReturn$VIP=VIP
  modelReturn$VIPPerRep=VIPRep
  # Average nVar over repetitions
  if (method=='PLS') {
    # Average nComp over repetitions
    nComp=mean(nCompRep)
    modelReturn$nComp=nComp
  }
  modelReturn$VAL$metric=fitness
  modelReturn$VAL$VAL=VAL
  if (modReturn) modelReturn$outModels=outMods
  modelReturn$yPredPerRep=yPred
  if (method=='PLS') modelReturn$nCompPerRep=nCompRep
  modelReturn$inData=InData
  ## Build overall "Fit" method for calculating R2 and visualisations
  if (method=='PLS'){
    if (DA) plsFit=plsda(X,Y,ncomp=round(nComp)) else 
      plsFit=pls(X,Y,ncomp=round(nComp),mode=methParam$mode)
    if (length(plsFit$nzv$Position)>0) removeVar=rownames(plsFit$nzv$Metrics) else removeVar=NA
    incVar=colnames(X)[!colnames(X)%in%removeVar]
    yFit=predict(plsFit,newdata=subset(X,select=incVar))$predict[,,round(nComp)]  # 
    modelReturn$Fit=list(yFit=yFit,plsFit=plsFit)
    # Prediction of newdata
  } else {
    rfFit=randomForest(X,Y)
    if (DA) {
      yFit=rfFit$votes
    } else {
      yFit=rfFit$predicted
    }
    modelReturn$Fit=list(yFit=yFit,rfFit=rfFit)
  }
  # Calculate fit statistics
  if (!DA) {
    TSS=sum((Y-mean(Y))^2)
    RSS=sum((Y-yFit)^2)
    PRESS=sum((Y-yPredAve)^2)
    R2=1-(RSS/TSS)
    Q2=1-(PRESS/TSS)
    modelReturn$fitMetric=data.frame(R2=R2,Q2=Q2)
  }
  # Stop timer
  end.time=proc.time()[3]
  modelReturn$calcMins=(end.time-start.time)/60
  cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
  class(modelReturn)=c('rdCVObject',method,ifelse(DA,'Classification',ifelse(ML,'Multilevel','Regression')))
  return(modelReturn)
}
