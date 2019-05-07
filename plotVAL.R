#' Plot validation metric 
#'
#' Produces a plot of validation metric vs number of variables in model (inner segment)
#' @param MVObj An object of class `MVObject`
#'
#' @return A plot
#' @export
plotVAL=function(MVObj) {
  VAL=MVObj$VAL$VAL
  metric=MVObj$VAL$metric
  count=as.numeric(colnames(VAL))
  nRep=dim(VAL)[3]
  plot(count,count,ylim=range(VAL),xlim=range(count),log='x',type='n',bty='l',ylab=metric,xlab='Number of variables (log scale)')
  for (r in 1:nRep) {
    matlines(count,t(VAL[,,r]),type='l',lty=1,col='lightgrey')
  }
  for (r in 1:nRep) {
    lines(count,colMeans(VAL[,,r]),col='darkgrey')
  }
  lines(count,apply(VAL,2,mean),col='black')
  for (i in 1:3) {
    abline(v=MVObj$nVar[i],lty=i,col=i+1,lwd=1.5)
  }
  legend('topleft',legend=c('Validation segments','Repetitions','Overall'),lty=1,col=c('lightgrey','darkgrey','black'),bty='n')
  legend('topright',legend=c("'Min (Minimal-optimal)","'Mid'","'Max' (All-relevant)"),lty=1:3,col=2:4,bty='n')
}
