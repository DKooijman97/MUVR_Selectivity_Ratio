#' Plot stability of selected variables and prediction fitness as a function of number of repetitions 
#'
#' @param MVObject MUVR object
#' @param model 'min' (default), 'mid' or 'max'
#' @param VAll Option of specifying which variables (i.e. names) to consider as reference set. Defaults to variables selected from the `model` of the `MVObject`
#' @param nVarLim Option of specifying upper limit for number of variables
#' @param missLim Option of specifying upper limit for number of misclassifications
#'
#' @return Plot of number of variables, proportion of variables overlapping with reference and prediction accuracy (Q2 for regression; MISS otherwise) as a function of number of repetitions. 
#' @export
plotStability=function(MVObject,model='min',VAll,nVarLim,missLim) {
  regr=any(class(MVObject)=='Regression')
  DA=MVObject$inData$DA
  ML=MVObject$inData$ML
  Y=MVObject$inData$Y
  nModel=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nVar=round(MVObject$nVar[nModel])
  if(missing(VAll)) VAll=names(sort(MVObject$VIP[,nModel])[1:nVar]) # Final selection of variables
  nRep=MVObject$inData$nRep
  nVRep=VARep=missRep=r2Rep=q2Rep=nV=VA=miss=r2=q2=numeric(nRep)
  for (i in 1:nRep) {
    nVRep[i]=MVObject$nVarPerRep[[nModel]][i]
    nV[i]=round(mean(MVObject$nVarPerRep[[nModel]][1:i]))
    VARep[i]=sum(names(sort(MVObject$VIPPerRep[[nModel]][,i])[1:nVRep[i]])%in%VAll)
    VA[i]=sum(names(sort(rowMeans(MVObject$VIPPerRep[[nModel]][,1:i,drop=F]))[1:nV[i]])%in%VAll)
    if(DA) {
      predsRep=MVObject$yPredPerRep[[nModel]][,,i,drop=F]
      missRep[i]=sum(levels(Y)[apply(predsRep,1,which.max)]!=Y)
      preds=MVObject$yPredPerRep[[nModel]][,,1:i]
      preds=apply(preds,c(1,2),mean)
      miss[i]=sum(levels(Y)[apply(preds,1,which.max)]!=Y)
    } else {
      predsRep=MVObject$yPredPerRep[[nModel]][,i,drop=F]
      PRESS=sum((Y-predsRep)^2)
      TSS=sum((Y-mean(Y))^2)
      q2Rep[i]=1-(PRESS/TSS)
      preds=MVObject$yPredPerRep[[nModel]][,1:i,drop=F]
      preds=rowMeans(preds)
      PRESS=sum((Y-preds)^2)
      TSS=sum((Y-mean(Y))^2)
      q2[i]=1-(PRESS/TSS)
    }
    if(ML) {
      class=ifelse(preds<0,-1,1)
      miss[i]=sum(class!=Y)
    }
  }
  VARep=VARep/length(VAll)
  VA=VA/length(VAll)
  if(missing(nVarLim)) {
    pot=10^floor(log10(max(nV)))
    nVarLim=ceiling(max(c(nV,nVRep))/pot)*pot
  }
  nPlot=ifelse(ML,4,3)
  par(mfrow=c(nPlot,1))
  par(mar=c(3,4,0,0)+.5)
  plot(nVRep,ylim=c(0,nVarLim),type='l',xlab='',ylab='Number of selected variables',col='grey',bty='l')
  lines(nV)
  legend('bottomright',c('Per repetition','Cumulative'),col=c('grey','black'),lty=1,bty='n')
  plot(VARep,type='l',ylim=c(0,1),col='pink',xlab='',ylab='Proportion of selected variables',bty='l')
  lines(VA,col='red')
  legend('bottomright',c('Per repetition','Cumulative'),col=c('pink','red'),lty=1,bty='n')
  if(DA | ML) {
    if(missing(missLim)) missLim=length(Y)
    plot(missRep,ylim=c(0,missLim),type='l',col='lightblue',xlab='',ylab='Number of misclassifications',bty='l')
    lines(miss,col='blue')
    legend('bottomright',c('Per repetition','Cumulative'),col=c('lightblue','blue'),lty=1,bty='n')
  }
  if(regr | ML) {
    plot(q2Rep,ylim=c(0,1),type='l',col='lightgreen',xlab='',ylab='Q2',bty='l')
    lines(q2,col='darkgreen')
    legend('bottomright',c('Per repetition','Cumulative'),col=c('lightgreen','darkgreen'),lty=1,bty='n')
  }
  mtext(text = 'Number of repetitions',side = 1,line = 2.3,cex=par()$cex)
  par(mfrow=c(1,1))
}
