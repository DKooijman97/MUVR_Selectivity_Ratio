#' MUVR: "Multivariate modelling with Unbiased Variable selection"
#' 
#' Repeated double cross validation with tuning of variables in the inner loop.
#'
#' @param X Predictor variables. NB: Variables (columns) must have names/unique identifiers. NAs not allowed in data. For multilevel, only the positive half of the difference matrix is specified.
#' @param Y Response vector (Dependent variable). For classification, a factor (or character) variable should be used. For multilevel, Y is calculated automatically.
#' @param ID Subject identifier (for sampling by subject; Assumption of independence if not specified)
#' @param scale If TRUE, the predictor variable matrix is scaled to unit variance for PLS modelling.
#' @param nRep Number of repetitions of double CV. (Defaults to 5)
#' @param nOuter Number of outer CV loop segments. (Defaults to 6)
#' @param nInner Number of inner CV loop segments. (Defaults to nOuter-1)
#' @param varRatio Ratio of variables to include in subsequent inner loop iteration. (Defaults to 0.75)
#' @param DA Boolean for Classification (discriminant analysis) (Defaults do FALSE, i.e. regression). PLS is limited to two-class problems (see `Y` above).
#' @param fitness Fitness function for model tuning (choose either 'AUROC' or 'MISS' (default) for classification; or 'RMSEP' (default) for regression.)
#' @param method Multivariate method. Supports 'PLS' and 'RF' (default)
#' @param methParam List with parameter settings for specified MV method (see function code for details)
#' @param ML Boolean for multilevel analysis (defaults to FALSE)
#' @param modReturn Boolean for returning outer segment models (defaults to FALSE)
#' @param nCompMax Option to choose max number of PLS components (default is 5)
#' @param logg Boolean for whether to sink model progressions to `log.txt`
#' @param parallel Boolean for whether to perform `foreach` parallel processing (Requires a registered parallel backend; Defaults to `TRUE`)
#'
#' @return A MUVR object
#' @export
MUVR=function(X,Y,ID,scale=TRUE,nRep=5,nOuter=6,nInner,varRatio=0.75,DA=FALSE,fitness=c('AUROC','MISS','BER','RMSEP'),method=c('PLS','RF'),nCompMax,methParam,ML=FALSE,modReturn=FALSE,logg=FALSE,parallel=TRUE){
  library(pROC)
  library(foreach)
  
  print("start MUVR")
  if (parallel) "%doVersion%"=get("%dopar%") else "%doVersion%"=get("%do%")
  if (missing(method)) method='RF'
  # Initialise modelReturn with function call
  modelReturn=list(call=match.call())
  # Start timer
  start.time=proc.time()[3]
  # Check indata
  if (length(dim(X))!=2) stop('\nWrong format of X matrix.\n')
  if (is.null(colnames(X))) stop('\nNo column names in X matrix.\n')
  X=as.matrix(X)
  if (method=='PLS') {
    nzv=MUVR::nearZeroVar(X)
    if (length(nzv$Position)>0) {
      modelReturn$nzv=colnames(X)[nzv$Position]
      X=X[,-nzv$Position]
      cat('\n',length(nzv$Position),'variables with near zero variance detected -> removed from X and stored under $nzv')
    }
  }
  nSamp=nrow(X)
  nVar=nVar0=ncol(X)
  if (missing(ID)) {
    cat('\nMissing ID -> Assume all unique (i.e. sample independence)')
    ID=1:nSamp
  }
  if (missing(nInner)) nInner=nOuter-1
  if (method=='RF') library(randomForest)
  if (missing(methParam)) {
    if (method=='PLS') {
      methParam=list(compMax=ifelse(nVar<5,nVar,5))
    } else {
      methParam=list(ntreeIn=150,ntreeOut=300,mtryMaxIn=150)
    }
    methParam$robust=0.05
  }
  if (!missing(nCompMax)) methParam$compMax=nCompMax
  if (ML) {
    X=rbind(X,-X)
    if (missing(Y)) Y=rep(-1,nSamp)
    Y=c(Y,-Y)
    nSamp=2*nSamp
    ID=c(ID,ID)
    DA=FALSE
    fitness='MISS'
    cat('\nMultilevel -> Regression on (-1,1) & fitness=MISS')
  }
  if (any(is.na(X)) | any(is.na(Y))) stop('\nNo missing values allowed in X or Y data.\n')
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
      fitness='MISS'
      cat('\nMissing fitness -> MISS')
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
  InData=list(X=X,Y=Y,ID=ID,scale=scale,nRep=nRep,nOuter=nOuter,nInner=nInner,varRatio=varRatio,DA=DA,fitness=fitness,method=method,methParam=methParam,ML=ML,parallel=parallel)
  ## Sort sampling based in subjects and not index
  unik=!duplicated(ID)  # boolean of unique IDs
  unikID=ID[unik]  
  if (DA) {
    if(nOuter>min(table(Y))) warning('\nnOuter is larger than your smallest group size. Consider lowering your nOuter to min(table(Y)).',call.=TRUE)
    unikY=Y[unik]  # Counterintuitive, but needed for groupings by Ynames
    Ynames=sort(unique(Y))  # Find groups
    groups=length(Ynames) # Number of groups
    groupID=list()  # Allocate list for indices of groups
    for (g in 1:groups) { 
      groupID[[g]]=unikID[unikY==Ynames[g]]  # Find indices per group
    }
    yPredMin=yPredMid=yPredMax=array(dim=c(length(Y),length(levels(Y)),nRep),dimnames=list(ID,levels(Y),paste('Rep',1:nRep,sep='')))
    yPredMinR=yPredMidR=yPredMaxR=matrix(nrow=length(Y),ncol=length(levels(Y)),dimnames=list(ID,levels(Y)))
  } else {
    yPredMin=yPredMid=yPredMax=matrix(nrow=length(Y),ncol=nRep,dimnames=list(ID,paste('Rep',1:nRep,sep='')))
    yPredMinR=yPredMidR=yPredMaxR=numeric(length(Y))
  }
  # Allocate response vectors and matrices for var's, nComp and VIP ranks over repetitions
  varRepMin=varRepMid=varRepMax=nCompRepMin=nCompRepMid=nCompRepMax=missRep=numeric(nRep)
  names(varRepMin)=names(varRepMid)=names(varRepMax)=names(nCompRepMin)=names(nCompRepMid)=names(nCompRepMax)=names(missRep)=paste(rep('rep',nRep),1:nRep,sep='')
  nCompSegMin=nCompSegMid=nCompSegMax=matrix(nrow=nRep,ncol=nOuter,dimnames=list(paste('repetition',1:nRep,sep=''),paste('segment',1:nOuter,sep='')))
  VIPRepMin=VIPRepMid=VIPRepMax=matrix(data=nVar0,nrow=nVar0,ncol=nRep,dimnames=list(colnames(X),paste(rep('rep',nRep),1:nRep,sep='')))
  var=numeric()
  cnt=0
  while (nVar>1) {  
    cnt=cnt+1
    var=c(var,nVar)
    nVar=floor(varRatio*nVar)
  }
  VAL=array(dim=c(nOuter,cnt,nRep),dimnames=list(paste('outSeg',1:nOuter,paste=''),var,paste(rep('rep',nRep),1:nRep,sep='')))
  ## Choose package/core algorithm according to chosen method
  packs=c('pROC')
  if(method=='RF') packs=c(packs,'randomForest')
  exports='vectSamp'
  ## Start repetitions
  # reps=list()
  # for (r in 1:nRep){
  reps=foreach(r=1:nRep, .packages=packs, .export=exports) %doVersion% {
    # r=1
    # r=r+1
    if (logg) sink('log.txt',append=TRUE)
    if (modReturn) outMod=list()
    cat('\n','   Repetition ',r,' of ',nRep,':',sep='')
    if (DA & identical(unikID,ID)) {
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
    varOutMin=varOutMid=varOutMax=nCompOutMin=nCompOutMid=nCompOutMax=numeric(nOuter)
    names(varOutMin)=names(varOutMid)=names(varOutMax)=names(nCompOutMin)=names(nCompOutMid)=names(nCompOutMax)=paste(rep('outSeg',nOuter),1:nOuter,sep='')
    VIPOutMin=VIPOutMid=VIPOutMax=matrix(data=nVar0,nrow=nVar0,ncol=nOuter,dimnames=list(colnames(X),paste(rep('outSeg',nOuter),1:nOuter,sep='')))
    VALRep=matrix(nrow=nOuter,ncol=cnt)
    ## Perform outer loop segments -> one "majority vote" MV model per segment
    for (i in 1:nOuter) {   
      # i=1
      # i=i+1
      cat('\n Segment ',i,' (variables):',sep='') # Counter
      ## Draw out test set
      testID=allTest[[i]] # Draw out segment = holdout set BASED ON UNIQUE ID
      testIndex=ID%in%testID # Boolean for samples included
      xTest=X[testIndex,]
      yTest=Y[testIndex]
      inID=unikID[!unikID%in%testID]  # IDs not in test set
      if (DA & identical(unikID,ID)) inY=unikY[!unikID%in%testID]  # Counterintuitive, but needed for grouping by Ynames
      ## Allocate variables for later use
      missIn=berIn=aucIn=rmsepIn=PRESSIn=nCompIn=matrix(nrow=nInner,ncol=cnt,dimnames=list(paste(rep('inSeg',nInner),1:nInner,sep=''),var))
      VIPInner=array(data=nVar0,dim=c(nVar0,cnt,nInner),dimnames=list(colnames(X),var,paste(rep('inSeg',nInner),1:nInner,sep='')))
      # Set variables
      incVar=colnames(X)
      ## Perform steps with successively fewer variables
      for (count in 1:cnt) {  # Build models with successively fewer variables. Quality metric = number of missclassifications for Validation set
        # count=1 # for testing
        # count=count+1
        nVar=var[count]
        cat(nVar)
        if (method=='PLS') comp=min(c(nVar,methParam$compMax))
        if (method=='RF') {
          mtryIn=ifelse(DA,min(c(methParam$mtryMaxIn,floor(sqrt(nVar)))),min(c(methParam$mtryMaxIn,floor(nVar/3))))
          mtryIn=max(c(2,mtryIn))
        }
        if (DA & identical(unikID,ID)) {
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
          xVal=subset(xVal,select=incVar)
          yVal=Y[valIndex]
          trainID=inID[!inID%in%valID]
          trainIndex=ID%in%trainID # Define Training segment
          xTrain=X[trainIndex,]
          xTrain=subset(xTrain,select=incVar)
          yTrain=Y[trainIndex]
          # sum(trainIndex,valIndex,testIndex)
          # trainIndex|valIndex|testIndex
          ## Make inner model
          if (method=='PLS') {
            source("plsInner.R")
            source("vip.R")
            inMod=plsInner(xTrain,yTrain,xVal,yVal,DA,fitness,comp,scale=scale)
            nCompIn[j,count]=inMod$nComp
          } else {
            inMod=MUVR::rfInner(xTrain,yTrain,xVal,yVal,DA,fitness,ntree=methParam$ntreeIn,mtry=mtryIn)
          }
          # Store fitness metric
          if (fitness=='MISS') {
            missIn[j,count]=inMod$miss
          } else if (fitness=='BER') {
            berIn[j,count]=inMod$ber
          } else if (fitness=='AUROC') {
            aucIn[j,count]=inMod$auc
          } else {
            rmsepIn[j,count]=inMod$rmsep
            PRESSIn[j,count]=(inMod$rmsep^2)*length(yVal)
          }
          # Store VIs
          VIPInner[match(names(inMod$vi),rownames(VIPInner)),count,j]=inMod$vi
        }
        ## Average inner VIP ranks before variable elimination
        VIPInAve=apply(VIPInner[,count,],1,mean)
        if (count<cnt) {
          incVar=names(VIPInAve[order(VIPInAve)])[1:var[count+1]]
        }
      }
      if (fitness=='AUROC') {
        fitRank=colMeans(-aucIn)
        VALRep[i,]=colMeans(aucIn)
      } else if (fitness=='MISS') {
        fitRank=VALRep[i,]=colSums(missIn)
      } else if (fitness=='BER') {
        fitRank=VALRep[i,]=colMeans(berIn)
      }else {
        fitRank=colMeans(rmsepIn)
        VALRep[i,]=sqrt(colSums(PRESSIn)/sum(!testIndex))
      }
      fitRank=(fitRank-min(fitRank))/abs(diff(range(fitRank))) # Rescale fitRank to range 0-1
      if(all(is.nan(fitRank))) fitRank=rep(0,cnt) # If all VAL have same value -> reset all fitRank to 0
      minIndex=max(which(fitRank<=methParam$robust))
      maxIndex=min(which(fitRank<=methParam$robust))
      # Per outer segment: Average inner loop variables, nComp and VIP ranks 
      varOutMin[i]=var[minIndex]
      varOutMax[i]=var[maxIndex]
      varOutMid[i]=round(exp(mean(log(c(var[minIndex],var[maxIndex])))))
      midIndex=which.min(abs(var-varOutMid[i]))
      if (method=='PLS') {
        nCompOutMin[i]=round(mean(nCompIn[,minIndex]))
        nCompOutMid[i]=round(mean(nCompIn[,midIndex]))
        nCompOutMax[i]=round(mean(nCompIn[,maxIndex]))
      }
      VIPOutMin[,i]=apply(VIPInner[,minIndex,],1,mean)
      VIPOutMid[,i]=apply(VIPInner[,midIndex,],1,mean)
      VIPOutMax[,i]=apply(VIPInner[,maxIndex,],1,mean)
      # Build outer model for min and max nComp and predict YTEST
      xIn=X[!testIndex,] # Perform Validation on all samples except holdout set
      yIn=Y[!testIndex]
      incVarMin=rownames(VIPOutMin)[rank(VIPOutMin[,i])<=varOutMin[i]]
      incVarMid=rownames(VIPOutMid)[rank(VIPOutMid[,i])<=varOutMid[i]]
      incVarMax=rownames(VIPOutMax)[rank(VIPOutMax[,i])<=varOutMax[i]]
      if (method=='PLS'){
        # Min model
        if (DA) plsOutMin=MUVR::plsda(subset(xIn,select=incVarMin),yIn,ncomp=nCompOutMin[i],near.zero.var=TRUE,scale=scale) else 
          plsOutMin=MUVR::pls(subset(xIn,select=incVarMin),yIn,ncomp=nCompOutMin[i],near.zero.var=TRUE,scale=scale)
        xTestMin=subset(xTest,select=incVarMin)
        yPredMinR[testIndex]=predict(plsOutMin,newdata=xTestMin,scale=scale)$predict[,,nCompOutMin[i]]  # 
        # Mid model
        if (DA) plsOutMid=MUVR::plsda(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale) else 
          plsOutMid=MUVR::pls(subset(xIn,select=incVarMid),yIn,ncomp=nCompOutMid[i],near.zero.var=TRUE,scale=scale)
        xTestMid=subset(xTest,select=incVarMid)
        yPredMidR[testIndex]=predict(plsOutMid,newdata=xTestMid,scale=scale)$predict[,,nCompOutMid[i]]  # 	
        # Max model
        if (DA) plsOutMax=MUVR::plsda(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale) else 
          plsOutMax=MUVR::pls(subset(xIn,select=incVarMax),yIn,ncomp=nCompOutMax[i],near.zero.var=TRUE,scale=scale)
        xTestMax=subset(xTest,select=incVarMax)
        yPredMaxR[testIndex]=predict(plsOutMax,newdata=xTestMax,scale=scale)$predict[,,nCompOutMax[i]]  # 
        if (modReturn) {
          outMod[[i]]=list(plsOutMin,plsOutMid,plsOutMax)
        }
      } else {
        rfOutMin=randomForest(x=subset(xIn,select=incVarMin),y=yIn,xtest=subset(xTest,select=incVarMin),ytest=yTest,keep.forest=TRUE)
        if (DA) {
          yPredMinR[testIndex,]=rfOutMin$test$votes
        } else {
          yPredMinR[testIndex]=rfOutMin$test$predicted
        }
        rfOutMid=randomForest(x=subset(xIn,select=incVarMid),y=yIn,xtest=subset(xTest,select=incVarMid),ytest=yTest,keep.forest=TRUE)
        if (DA) {
          yPredMidR[testIndex,]=rfOutMid$test$votes
        } else {
          yPredMidR[testIndex]=rfOutMid$test$predicted
        }
        rfOutMax=randomForest(x=subset(xIn,select=incVarMax),y=yIn,xtest=subset(xTest,select=incVarMax),ytest=yTest,keep.forest=TRUE)
        if (DA) {
          yPredMaxR[testIndex,]=rfOutMax$test$votes
        } else {
          yPredMaxR[testIndex]=rfOutMax$test$predicted
        }
        if (modReturn) {
          outMod[[i]]=list(rfOutMin,rfOutMid,rfOutMax)
        }
      }
    }
    # Per repetition: Average outer loop predictions, VIP ranks and nComp for PLS
    parReturn=list(yPredMin=yPredMinR,yPredMid=yPredMidR,yPredMax=yPredMaxR)
    parReturn$VIPRepMin=apply(VIPOutMin,1,mean)
    parReturn$VIPRepMid=apply(VIPOutMid,1,mean)
    parReturn$VIPRepMax=apply(VIPOutMax,1,mean)
    if (method=='PLS'){
      parReturn$nCompRepMin=round(mean(nCompOutMin))
      parReturn$nCompRepMid=round(mean(nCompOutMid))
      parReturn$nCompRepMax=round(mean(nCompOutMax))
      parReturn$nCompSegMin=nCompOutMin
      parReturn$nCompSegMid=nCompOutMid
      parReturn$nCompSegMax=nCompOutMax
    }
    # Validation curves
    parReturn$VAL=VALRep 
    # Calculate nVar per repetition
    fitRankRep=colSums(VALRep)
    if(fitness=='AUROC') fitRankRep=-fitRankRep
    fitRankRep=(fitRankRep-min(fitRankRep))/abs(diff(range(fitRankRep)))
    if(all(is.nan(fitRankRep))) fitRankRep=rep(0,cnt) # If all VAL have same value -> reset all fitRankRep to 0
    minIndex=max(which(fitRankRep<=methParam$robust))
    maxIndex=min(which(fitRankRep<=methParam$robust))
    parReturn$varRepMin=var[minIndex]
    parReturn$varRepMid=round(exp(mean(log(c(var[minIndex],var[maxIndex])))))
    parReturn$varRepMax=var[maxIndex]
    # Return underlying models
    if (modReturn) parReturn$outModel=outMod
    if (logg) sink()
    return(parReturn)
    # reps[[r]]=parReturn
  }
  if (modReturn) outMods=list()
  for (r in 1:nRep) {
    if (DA) yPredMin[,,r]=reps[[r]]$yPredMin else
      yPredMin[,r]=reps[[r]]$yPredMin
    if (DA) yPredMid[,,r]=reps[[r]]$yPredMid else
      yPredMid[,r]=reps[[r]]$yPredMid
    if (DA) yPredMax[,,r]=reps[[r]]$yPredMax else
      yPredMax[,r]=reps[[r]]$yPredMax
    varRepMin[r]=reps[[r]]$varRepMin
    varRepMid[r]=reps[[r]]$varRepMid
    varRepMax[r]=reps[[r]]$varRepMax
    VIPRepMin[,r]=reps[[r]]$VIPRepMin
    VIPRepMid[,r]=reps[[r]]$VIPRepMid
    VIPRepMax[,r]=reps[[r]]$VIPRepMax
    if (method=='PLS') {
      nCompRepMin[r]=reps[[r]]$nCompRepMin
      nCompRepMid[r]=reps[[r]]$nCompRepMid
      nCompRepMax[r]=reps[[r]]$nCompRepMax
      nCompSegMin[r,]=reps[[r]]$nCompSegMin
      nCompSegMid[r,]=reps[[r]]$nCompSegMid
      nCompSegMax[r,]=reps[[r]]$nCompSegMax
    }
    VAL[,,r]=reps[[r]]$VAL
    if (modReturn) outMods=c(outMods,reps[[r]]$outModel)
  }
  # Average predictions
  if (DA) {
    yPred=list()
    yPred[['min']]=apply(yPredMin,c(1,2),mean)
    yPred[['mid']]=apply(yPredMid,c(1,2),mean)
    yPred[['max']]=apply(yPredMax,c(1,2),mean)
  } else {
    yPred=apply(yPredMin,1,mean)[1:nrow(X)]
    yPred=cbind(yPred,apply(yPredMid,1,mean)[1:nrow(X)])
    yPred=cbind(yPred,apply(yPredMax,1,mean)[1:nrow(X)])
    colnames(yPred)=c('min','mid','max')
    rownames(yPred)=paste(1:nSamp,ID,sep='_ID')
  }
  modelReturn$yPred=yPred
  if (DA) {
    auc=matrix(nrow=3,ncol=length(levels(Y)),dimnames=list(c('min','mid','max'),levels(Y)))
    for (cl in 1:length(levels(Y))) {
      auc[1,cl]=roc(Y==(levels(Y)[cl]),yPred[['min']][,cl])$auc
      auc[2,cl]=roc(Y==(levels(Y)[cl]),yPred[['mid']][,cl])$auc
      auc[3,cl]=roc(Y==(levels(Y)[cl]),yPred[['max']][,cl])$auc
    }
    # Classify predictions
    miss=numeric(3)
    yClass=data.frame(Y)
    for (mo in 1:3) {
      classPred=factor(apply(yPred[[mo]],1,function(x) levels(Y)[which.max(x)]),levels=levels(Y))
      miss[mo]=sum(classPred!=Y)
      yClass[,mo]=classPred
    }
    names(miss)=colnames(yClass)=c('min','mid','max')
    rownames(yClass)=paste(1:nSamp,ID,sep='_ID')
    # Report
    modelReturn$yClass=yClass
    modelReturn$miss=miss
    modelReturn$auc=auc
  } else if (ML) {
    modelReturn$yClass=apply(yPred,2,function(x) ifelse(x>0,1,-1))
    modelReturn$miss=apply(modelReturn$yClass,2,function(x) sum(x!=Y))
    modelReturn$auc=apply(yPred,2,function(x) roc(Y,x)$auc)
    colnames(modelReturn$yClass)=names(modelReturn$miss)=names(modelReturn$auc)=c('min','mid','max')
    rownames(modelReturn$yClass)=paste(1:nSamp,ID,sep='_ID')
  }
  # Average VIP ranks over repetitions
  VIP=apply(VIPRepMin,1,mean)
  VIP=cbind(VIP,apply(VIPRepMid,1,mean))
  VIP=cbind(VIP,apply(VIPRepMax,1,mean))
  colnames(VIP)=c('min','mid','max')
  modelReturn$VIP=VIP
  modelReturn$VIPPerRep=list(minModel=VIPRepMin,midModel=VIPRepMid,maxModel=VIPRepMax)
  # Calculate overall nVar
  fitRankAll=apply(VAL,2,mean)
  if(fitness=='AUROC') fitRankAll=-fitRankAll
  fitRankAll=(fitRankAll-min(fitRankAll))/abs(diff(range(fitRankAll))) # rescale to 0-1 range
  if(all(is.nan(fitRankAll))) fitRankAll=rep(0,cnt) # If all VAL have same value -> reset all fitRankAll to 0
  minIndex=max(which(fitRankAll<=methParam$robust))
  maxIndex=min(which(fitRankAll<=methParam$robust))
  nVar=c(var[minIndex],round(exp(mean(log(c(var[minIndex],var[maxIndex]))))),var[maxIndex])
  names(nVar)=c('min','mid','max')
  modelReturn$nVar=nVar
  if (method=='PLS') {
    # Average nComp over repetitions
    nComp=c(round(mean(nCompRepMin)),round(mean(nCompRepMid)),round(mean(nCompRepMax)))
    names(nComp)=c('min','mid','max')
    modelReturn$nComp=nComp
  }
  modelReturn$VAL$metric=fitness
  modelReturn$VAL$VAL=VAL
  if (modReturn) modelReturn$outModels=outMods
  modelReturn$yPredPerRep=list(minModel=yPredMin,midModel=yPredMid,maxModel=yPredMax)
  modelReturn$nVarPerRep=list(minModel=varRepMin,midModel=varRepMid,maxModel=varRepMax)
  if (method=='PLS') {
    modelReturn$nCompPerRep=list(minModel=nCompRepMin,midModel=nCompRepMid,maxModel=nCompRepMax)
    modelReturn$nCompPerSeg=list(minModel=nCompSegMin,midModel=nCompSegMid,maxModel=nCompSegMax)
  }
  modelReturn$inData=InData
  ## Build overall "Fit" method for calculating R2 and visualisations
  incVarMin=names(VIP[rank(VIP[,1])<=round(nVar[1]),1])
  incVarMid=names(VIP[rank(VIP[,2])<=round(nVar[2]),2])
  incVarMax=names(VIP[rank(VIP[,3])<=round(nVar[3]),3])
  if (method=='PLS'){
    # Min model
    if (DA) plsFitMin=MUVR::plsda(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale) else 
      plsFitMin=MUVR::pls(subset(X,select=incVarMin),Y,ncomp=round(nComp[1]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMin$nzv$Position)>0) incVarMin=incVarMin[!incVarMin%in%rownames(plsFitMin$nzv$Metrics)]
    yFitMin=predict(plsFitMin,newdata=subset(X,select=incVarMin),scale=scale)$predict[,,nComp[1]]  # 
    # Mid model
    if (DA) plsFitMid=MUVR::plsda(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale) else 
      plsFitMid=MUVR::pls(subset(X,select=incVarMid),Y,ncomp=round(nComp[2]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMid$nzv$Position)>0) incVarMid=incVarMid[!incVarMid%in%rownames(plsFitMid$nzv$Metrics)]
    yFitMid=predict(plsFitMid,newdata=subset(X,select=incVarMid),scale=scale)$predict[,,nComp[2]]  # 
    # Max model
    if (DA) plsFitMax=MUVR::plsda(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale) else 
      plsFitMax=MUVR::pls(subset(X,select=incVarMax),Y,ncomp=round(nComp[3]),near.zero.var=TRUE,scale=scale)
    if (length(plsFitMax$nzv$Position)>0) incVarMax=incVarMax[!incVarMax%in%rownames(plsFitMax$nzv$Metrics)]
    yFitMax=predict(plsFitMax,newdata=subset(X,select=incVarMax),scale=scale)$predict[,,nComp[3]]  # 
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,plsFitMin=plsFitMin,plsFitMid=plsFitMid,plsFitMax=plsFitMax)
  } else {
    rfFitMin=suppressWarnings(randomForest(subset(X,select=incVarMin),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMin=rfFitMin$votes
    } else {
      yFitMin=rfFitMin$predicted
    }
    rfFitMid=suppressWarnings(randomForest(subset(X,select=incVarMid),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMid=rfFitMid$votes
    } else {
      yFitMid=rfFitMid$predicted
    }
    rfFitMax=suppressWarnings(randomForest(subset(X,select=incVarMax),Y)) # suppress warnings for ML (regression against fewer than 5 unique values)
    if (DA) {
      yFitMax=rfFitMax$votes
    } else {
      yFitMax=rfFitMax$predicted
    }
    yFit=cbind(yFitMin,yFitMid,yFitMax)
    yRep=ncol(yFit)/3
    colnames(yFit)=rep(c('min','mid','max'),each=yRep)
    rownames(yFit)=ID
    modelReturn$Fit=list(yFit=yFit,rfFitMin=rfFitMin,rfFitMid=rfFitMid,rfFitMax=rfFitMax)
  }
  # Calculate fit statistics
  if (!DA) {
    TSS=sum((Y-mean(Y))^2)
    RSSMin=sum((Y-yFitMin)^2)
    RSSMid=sum((Y-yFitMid)^2)
    RSSMax=sum((Y-yFitMax)^2)
    PRESSMin=sum((Y-yPred[,1])^2)
    PRESSMid=sum((Y-yPred[,2])^2)
    PRESSMax=sum((Y-yPred[,3])^2)
    R2Min=1-(RSSMin/TSS)
    R2Mid=1-(RSSMid/TSS)
    R2Max=1-(RSSMax/TSS)
    Q2Min=1-(PRESSMin/TSS)
    Q2Mid=1-(PRESSMid/TSS)
    Q2Max=1-(PRESSMax/TSS)
    modelReturn$fitMetric=list(R2=c(R2Min,R2Mid,R2Max),Q2=c(Q2Min,Q2Mid,Q2Max))
  }
  # Stop timer
  end.time=proc.time()[3]
  modelReturn$calcMins=(end.time-start.time)/60
  cat('\n Elapsed time',(end.time-start.time)/60,'mins \n')
  class(modelReturn)=c('MVObject',method,ifelse(DA,'Classification',ifelse(ML,'Multilevel','Regression')))
  return(modelReturn)
}

