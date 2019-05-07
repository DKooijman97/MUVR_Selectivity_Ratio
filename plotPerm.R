#' Plot for comparison of actual model fitness vs permutation
#' 
#' Plots histogram of null hypothesis (permutation) distribution, actual model fitness and cumulative p-value.
#' Plot defaults to "greater than" or "smaller than" tests and cumulative probability in Student's t-distribution
#'
#' @param actual Actual model fitness (e.g. Q2, AUROC or number of misclassifications)
#' @param h0 Null hypothesis (permutation) distribution of similar metric as `actual`
#' @param xlab Label for x-axis (e.g. 'Q2', 'AUROC', or 'Misclassifications')
#' @param side Cumulative p either "greater" or "smaller" than H0 distribution (defaults to side of median(H0))
#' @param type Choice of Student's t-distribution of original ('t', default) or ranked ('non') data for non-parametric test
#' @param xlim Choice of user-specified x-limits (if default is not adequate)
#' @param ylim Choice of user-specified y-limits (if default is not adequate)
#' @param breaks Choice of user-specified histogram breaks (if default is not adequate)
#' @param main Choice of user-specified plot title
#' @param pos Choice of position of p-value label (if default is not adequate)
#'
#' @return Plot
#' @export
plotPerm=function(actual,h0,xlab=NULL,side=c('greater','smaller'),type=c('t','non'),xlim,ylim=NULL,breaks='Sturges',pos, main=NULL) {
  if(missing(side)) side=ifelse(actual<median(h0),'smaller','greater')
  if(missing(type)) type='t'
  if(missing(pos)) pos=ifelse(side=='smaller',4,2)
  pP=pPerm(actual,h0,side,type)
  if(missing(xlim)) {
    if(side=='smaller') xlim=c(0,max(h0)) else xlim=c(min(h0),1) 
  }
  (h=hist(h0,breaks,xlim=xlim,ylim=ylim,axes=F,xlab=xlab,freq=FALSE,main=main))
  h2=max(h$density)*.75
  axis(1,pos=0)
  if(side=='smaller') axis(2,pos=0,las=1) else axis(2,pos=h$breaks[1],las=1)
  lines(rep(actual,2),c(0,h2))
  text(actual,h2,pos=pos,labels=paste('p=',signif(pP,4),sep=''))
}
