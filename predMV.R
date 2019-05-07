#' Predict MV object
#'
#' At present, this function only supports predictions for PLS regression type problems
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param newdata New data for which to predict outcomes
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#'
#' @return A pdf with plots of results from multivariate predictions
#' @export
predMV=function(MVObj,newdata,model='mid') {
  if (!any(class(MVObj)=='MVObject')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  modNum=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  method=MVObj$inData$method
  nRep=MVObj$inData$nRep
  nOuter=MVObj$inData$nOuter
  if (method=='PLS') {
    nComps=MVObj$nCompPerSeg[[modNum]]
  } else {
    library(randomForest)
  }
  # par(mar=c(4,4,0,0)+.5)
  if (class(MVObj)[3]=='Regression') {
    yPredPerMod=matrix(ncol=length(MVObj$outModels),nrow=nrow(newdata),dimnames=list(paste('observation',1:nrow(newdata),sep=''),paste('model',1:length(MVObj$outModels),sep='')))
    n=0
    for(r in 1:nRep) {
      for(i in 1:nOuter) {
        n=n+1
        mod=MVObj$outModels[[n]][[modNum]]
        if (method=='PLS') {
          if (any(!colnames(mod$X)%in%colnames(newdata))) {
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {
            X=subset(newdata,select=colnames(mod$X))
            nComp=nComps[r,i]
            yPredPerMod[,n]=predict(mod,newdata=X)$predict[,,nComp]  # 
          }
        } else {
          if (any(!rownames(mod$importance)%in%colnames(newdata))) {
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {
            X=subset(newdata,select=mod$names$X)
            X=subset(newdata,select=rownames(mod$importance))
            yPredPerMod[,n]=predict(mod,newdata=X)  # 
          }
        }
      }
    }
    yPred=apply(yPredPerMod,1,mean)
    return(list(yPred=yPred,yPredPerMod=yPredPerMod))
  } else if (class(MVObj)[3]=='Classification') {
    yPredPerMod=array(dim=c(nrow(newdata),length(levels(MVObj$inData$Y)),length(MVObj$outModels)),dimnames=list(paste('observation',1:nrow(newdata),sep=''),levels(MVObj$inData$Y),paste('model',1:length(MVObj$outModels),sep='')))
    n=0
    for(r in 1:nRep) {
      for(i in 1:nOuter) {
        n=n+1
        mod=MVObj$outModels[[n]][[modNum]]
        if (method=='PLS') {
          cat('\nNot yet implemented')
          return(NULL)
        } else {
          if (any(!rownames(mod$importance)%in%colnames(newdata))) {
            cat('\nMismatch variable names in model',n,': Return NULL')
            return(NULL)
          } else {
            X=subset(newdata,select=rownames(mod$importance))
            yPredPerMod[,,n]=predict(mod,newdata=X,type='vote')  # 
          }
        }
      }
    }
    yPred=apply(yPredPerMod,c(1,2),mean)
    yClass=levels(MVObj$inData$Y)[apply(yPred,1,which.max)]
    names(yClass)=paste('observation',1:nrow(newdata),sep='')
    return(list(yClass=yClass,yPred=yPred,yPredPerMod=yPredPerMod))
  } else {
    cat('\nNot yet implemented')
  }
}
