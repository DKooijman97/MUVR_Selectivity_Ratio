#' Plot predictions
#'
#' At present, this function only supports predictions for PLS regression type problems
#'
#' @param Ytrue 
#' @param Ypreds 
#'
#' @export
plotPred=function(Ytrue,Ypreds) {
  par(mar=c(4,4,0,0)+.5)
  matplot(Ytrue,t(Ypreds),pch='.',col='grey')
  points(Ytrue,colMeans(Ypreds),pch='.',col='black',cex=2)
}
