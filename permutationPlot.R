#' Plot permutation analysis using actual model and permutations()
#'
#' This is asically a wrapper for the MUVR::plotPerm() function using model objects to make coding nicer and cleaner
#'
#' @param MVObj A MUVR object
#' @param permMatrix A permutation fitness outcome from MUVR::permutations()
#' @param model 'Min', 'Mid', or 'Max'
#' @param type 't' (default; for Student's t) or 'non' for "non-parametric" (i.e. rank) studen'ts
#' @param side 'smaller' for actual lower than H0 or 'greater' for actual larger than H0 (automatically selected if not specified)
#' @param pos which side of actual to put p-value on
#' @param xlab optional xlabel
#' @param xlim optional x-range
#' @param ylim otional y-range
#' @param breaks optional custom histogram breaks (defaults to 'sturges')
#' @param main optional plot title (or TRUE for autoname)
#'
#' @return A permutation plot
#' @export
#'
#' @examples
#' library(MUVR)
#' library(doParallel)
#' nCore=detectCores()-1
#' cl=makeCluster(nCore)
#' registerDoParallel(cl)
#' nRep=2*nCore
#' varRatio=.75
#' nOuter=6
#' nPerm=50
#' R12ML=MUVR(X=mlr12,ML=T,nRep=nRep,nOuter=nOuter,varRatio=varRatio,method='RF')
#' permR12=permutations(R12ML)
#' stopCluster(cl)
#' permutationPlot(R12ML,permR12)
permutationPlot=function(MVObj,permMatrix,model='Mid',type=c('t','non'),side=c('greater','smaller'),pos,xlab=NULL,xlim,ylim=NULL,breaks='Sturges',main=NULL) {
  nModel=ifelse(model=='Min',1,ifelse(model=='Mid',2,3))
  if (any(class(MVObj)=='Regression')) {
    actual=MVObj$fitMetric$Q2[nModel] 
    if (missing(xlab)) xlab='Q2'
  } else {
    actual=MVObj$miss[nModel]
    if (missing(xlab)) xlab='Misclassifications'
  }
  h0=permMatrix[,nModel]
  if(missing(side)) side=ifelse(actual<median(h0),'smaller','greater')
  if(missing(type)) type='t'
  if(missing(pos)) pos=ifelse(side=='smaller',4,2)
  if(missing(xlim)) {
    if(side=='smaller') xlim=c(0,max(h0)) else xlim=c(min(h0),1) 
  }
  if(isTRUE(main)) main=paste('Permutation analysis of',deparse(substitute(MVObj)))
  plotPerm(actual = actual, h0 = h0, type = type, pos = pos, side = side, xlab = xlab, xlim = xlim, ylim = ylim, breaks = breaks, main = main)
}
