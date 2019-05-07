#' PLS model in inner CV loop
#'
#' Called from Wrapper
#'
#' @param xTrain Training data (samples as rows; variables as columns)
#' @param yTrain Training response
#' @param xVal Validation data
#' @param yVal Validation response (for tuning)
#' @param DA Logical for discriminant analysis (classification)
#' @param fitness Fitness function ('MISS', 'AUROC' or 'RMSEP')
#' @param comp Max number of components to try
#' @param scale Whether or not to scale inData (X)
#'
#' @return An object containing:
#' @return (`miss`, `auc` or `rmsep`) A fitness metric
#' @return `nComp` Optimised number of components within range (1:comp)
#' @return `vi` variable importance rankings
#' @export
#'
plsInner=function(xTrain,yTrain,xVal,yVal,DA,fitness,comp,scale=TRUE) {
  cond=TRUE
  while(cond) {
    yValInner=tryCatch({
      if (DA) plsModIn=MUVR::plsda(xTrain,yTrain,ncomp=comp,near.zero.var=TRUE,scale=scale) else plsModIn=MUVR::pls(xTrain,yTrain,ncomp=comp,near.zero.var=TRUE,scale=scale)
      yValInner=predict(plsModIn,newdata=xVal,onlyPred=TRUE,scale=scale)$predict[,,]
    }, error=function(e) return('error'))
    if ((length(yValInner)==1 && yValInner=='error') | any(is.na(yValInner))) comp=comp-1 else cond=FALSE
    if (comp==0) cond=FALSE
  }
  returnIn=list()
  if (comp>0){
    if(!DA & !is.matrix(yValInner)) yValInner=as.matrix(yValInner)
    if (DA) {
      if (fitness=='MISS') {
        if(comp>1) classes=apply(yValInner,c(1,3),which.max) else classes=matrix(apply(yValInner,1,which.max),ncol=1)
        misClass=apply(classes,2,function(x) sum(x!=as.numeric(yVal)))
        returnIn$miss=min(misClass,na.rm=T)
        nComp=which.min(misClass)
      } else if (fitness=='BER') {
        if(comp>1) classes=apply(yValInner,c(1,3),which.max) else classes=matrix(apply(yValInner,1,which.max),ncol=1)
        BER=apply(classes,2,function(x) getBER(actual=as.numeric(yVal),predicted=x))
        returnIn$ber=min(BER,na.rm=T)
        nComp=which.min(BER)
      } else {
        auc=apply(yValInner[,1,],2,function(x) roc(yVal,x)$auc)
        returnIn$auc=max(auc,na.rm=T)
        nComp=which.max(auc)
      }
    } else {
      if (fitness=='MISS') {
        # cat(' miss',count)
        yClassInner=ifelse(yValInner>0,1,-1)
        misClass=apply(yClassInner,2,function(x) sum(x!=yVal))
        returnIn$miss=min(misClass,na.rm=T)
        nComp=which.min(misClass)
      }
      if (fitness=='AUROC') {
        # cat(' auc',count)
        auc=apply(yValInner,2,function(x) roc(yVal,x)$auc)
        returnIn$auc=max(auc,na.rm=T)
        nComp=which.max(auc)
      }
      if (fitness=='RMSEP') {
        # cat(' rmsep',count)
        rmsep=apply(yValInner,2,function(x) sqrt(sum((yVal-x)^2,na.rm=T)/(length(yVal)-sum(is.na(x)))))
        returnIn$rmsep=min(rmsep,na.rm=T)
        nComp=which.min(rmsep)
      }
    }
    #source("vip.R")
    returnIn$vi=rank(-vip(plsModIn)[,nComp])    #Call to VIP, plsMod is the object [, nComp] creates vector of 1 to n PC. DennisK 
  } else {
    nComp=0
    if (fitness=='MISS') returnIn$miss=length(yVal)   #Multiple fitness tests?  DennisK
    if (fitness=='AUROC') returnIn$auc=0
    if (fitness=='BER') returnIn$ber=1
    if (fitness=='RMSEP') returnIn$rmsep=1E10
    returnIn$vip=rep(1,ncol(xTrain))                  #Temporary save of VIPS in inner loop? DennisK 
  }
  returnIn$nComp=nComp
  return(returnIn)
}
