#' Make permutations with data and default settings from an actual MUVR object
#' 
#' This function will extract data and parameter settings from a MUVR object and run standard permutations.
#' This will fit a standard case of multivariate predictive modelling in either a regression, classification or multilevel case.
#' However, if an analysis has a complex sample dependency which requires constrained permutation of your response vector 
#' or if a variable pre-selection is performed for decreased computational burden, then permutaion loops should be constructed manually.
#' In those cases, View(permutations) can be a first start from which to build custom solutions for permutation analysis.
#' @param MVObj A MUVR obvject
#' @param nPerm number of permutations to run
#' @param nRep number of repetitions for each permutation (defaults to value of actual model)
#' @param nOuter number of outer validation segments for each permutation (defaults to value of actual model)
#' @param varRatio varRatio for each permutation (defaults to value of actual model)
#' @param parallel whether to run calculations using parallel processing - which requires registered backend (defaults to value of actual model)
#'
#' @return A permutation matrix with permuted fitness statistics (nrow=nPerm and ncol=3 for min/mid/max)
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

permutations=function(MVObj,nPerm=50,nRep,nOuter,varRatio,parallel) {
  name=deparse(substitute(MVObj))
  X=MVObj$inData$X
  Y=MVObj$inData$Y
  ID=MVObj$inData$ID
  DA=MVObj$inData$DA
  ML=MVObj$inData$ML
  scale=MVObj$inData$scale
  fitness=MVObj$inData$fitness
  method=MVObj$inData$method
  methParam=MVObj$inData$methParam
  if(missing(nRep)) nRep=MVObj$inData$nRep
  if(missing(nOuter)) {
    nOuter=MVObj$inData$nOuter
    nInner=MVObj$inData$nInner
  } else nInner=nOuter-1
  if(missing(varRatio)) varRatio=MVObj$inData$varRatio
  if(missing(parallel)) parallel=MVObj$inData$parallel
  if(ML) {
    nSamp=nrow(X)/2
    X=X[1:nSamp,]
    ID=ID[1:nSamp]
  }
  startTime=proc.time()[3]
  h0=matrix(ncol=3,nrow=nPerm)
  colnames(h0)=c('Min','Mid','Max')
  for (p in 1:nPerm) {
    cat('\n"',name,'" permutation ',p,' of ',nPerm,'\n',sep = '')
    if (ML) YPerm=sample(c(-1,1),size=nSamp,replace=TRUE) else YPerm=sample(Y)
    permMod=MUVR(X=X, Y=YPerm, ID=ID, scale=scale, DA=DA, ML=ML, nRep=nRep, nOuter=nOuter, nInner=nInner, varRatio=varRatio, fitness=fitness, method=method, methParam=methParam, parallel=parallel)
    if (any(class(MVObj)=='Regression')) h0[p,]=permMod$fitMetric$Q2 else h0[p,]=permMod$miss
    nowTime=proc.time()[3]
    timePerRep=(nowTime-startTime)/p
    timeLeft=(timePerRep*(nPerm-p))/60
    cat('\nEstimated time left:',timeLeft,'mins\n\n')
  }
  return(h0)
}
