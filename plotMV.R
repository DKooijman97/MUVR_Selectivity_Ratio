#' Plot MV object
#'
#' @param MVObj An `MVobject` obtained from the MVWrap function
#' @param model What type of model to plot ('min', 'mid' or 'max'). Defaults to 'mid'.
#' @param factCols An optional vector with colors for the factor levels (in the same order as the levels)
#' @param sampLabels Sample labels (optional; implemented for classification)
#' @param ylim Optional for imposing y-limits for regression and classification analysis
#'
#' @return A plot of results from multivariate predictions
#' @export
plotMV=function(MVObj,model='mid',factCols,sampLabels,ylim=NULL) {
  if (!any(class(MVObj)=='MVObject')) {
    cat('\nWrong object class: Return NULL')
    return(NULL)
  }
  modNum=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  Y=MVObj$inData$Y
  nSamp=length(Y)
  if(missing(sampLabels)) sampLabels=Y
  if(length(sampLabels)!=nSamp) {
    warning('Length of sampLabels not equal to number of samples in Y. \n   Autonumbers used instead.')
    sampLabels=1:nSamp
  }
  # par(mar=c(4,4,2,0)+.5)
  if (class(MVObj)[3]=='Regression') {
    YP=MVObj$yPred[,modNum]
    YPR=MVObj$yPredPerRep[[modNum]]
    if (is.null(ylim)) ylim <- range(YPR)
    matplot(Y,YPR,pch=20,xlab='Original Y',ylab='Predicted Y',col='grey',bty='l',cex=0.5, ylim=ylim)
    points(Y,YP,pch=20)
    reg=lm(YP~Y)
    abline(reg)
    legend('topleft',legend=c(paste('Model R2 =',signif(MVObj$fitMetric$R2[modNum],3)),paste('Model Q2 =',signif(MVObj$fitMetric$Q2[modNum],3))),bty='n')
  } else if (class(MVObj)[3]=='Classification') {
    YP=MVObj$yPred[[modNum]]
    YPR=MVObj$yPredPerRep[[modNum]]
    if (is.null(ylim)) ylim <- range(YPR)
    classes=1:length(levels(Y))
    if(missing(factCols)) factCols=classes+1
    if(length(factCols)!=length(classes)) {
      warning('Length of factCols not equal to number of levels in Y. \n   Autocolors used instead.')
      factCols=classes+1
    }
    classNudge=0.2*((classes-mean(classes))/(mean(classes)-1))
    plot(1:nSamp,Y,type='n',ylim=ylim,xlab='',ylab='Class prediction probability',xaxt='n')
    axis(1,at=1:length(Y),labels = sampLabels,las=3)
    for(cl in classes) {
      matpoints((1:nSamp)+classNudge[cl],YPR[,cl,],pch=20,col=factCols[cl],cex=0.5)
      points((1:nSamp)+classNudge[cl],YP[,cl],pch=20,col=factCols[cl])
    }
    for (li in 1:(nSamp+1)) {
      abline(v=li-.5,lty=3,col='grey')
    }
    yClass=MVObj$yClass[,modNum]
    whichWrong=which(yClass!=Y)
    wrongClass=as.numeric(Y[whichWrong])
    for (w in 1:length(wrongClass)) {
      points(whichWrong[w]+classNudge[wrongClass[w]],YP[whichWrong[w],wrongClass[w]],cex=2)
    }
    legend('topleft',legend=c(levels(Y),'misclassified'),pch=c(rep(16,length(classes)),1),col=c(factCols,1),cex=0.8,pt.cex=c(rep(0.5,length(classes)),2),bty='n')
  } else {
    YP=MVObj$yPred[,modNum]
    YPR=MVObj$yPredPerRep[[modNum]]
    matplot(YPR,1:nSamp,pch=20,col='grey',cex=0.5,ylim=c(nSamp,1),ylab='Sample number',xlab='Predicted Y')
    points(YP,1:nSamp,pch=20,col='black')
    abline(h=nSamp/2+0.5,lty=2)
    abline(v=0,lty=2)
  }
}
