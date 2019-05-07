#' Perform matrix pre-processing
#'
#' @param X Data matrix with samples in rows and variables in columns
#' @param offset Add offset to all data points (defaults to 0)
#' @param zeroOffset Add offset to zero data (defaults to 0)
#' @param trans Either 'log', 'sqrt' or 'none' (default is 'none')
#' @param center Either 'mean', 'none' or a numeric vector of length equal to the number of columns of X (defaults to 'none').
#' @param scale Either 'UV', 'Pareto', 'none' or a numeric vector of length equal to the number of columns of X (defaults to 'none').
#'
#' @return A pre-processed data matrix
#' @export
preProcess=function(X,offset=0,zeroOffset=0,trans='none',center='none',scale='none') {
  nVar=ncol(X)
  # Add offset
  X[X==0]=zeroOffset
  cat('Zero offset:',zeroOffset)
  X=X+offset
  cat('\nOffset:',offset)
  # Perform transformation
  trans=match.arg(trans,c('log','sqrt','none'))
  cat('\nTransformation:',trans)
  if(trans!='none') {
    if(trans=='log') {
      if(any(X<=0)) stop('no zero or negative values allowed when doing log transformation')
      X=apply(X,2,log) 
    }
    if(trans=='sqrt') {
      if(any(X<0)) stop('no negative values allowed when doing sqrt transformation')
      X=apply(X,2,sqrt) 
    }
  }
  # Perform centering and scaling
  if (length(center)!=nVar) {
    center=match.arg(center,c('mean','none'))
    cat('\nCenter:',center)
    if (center=='mean') center=TRUE else center=FALSE
  } else cat('\nCenter: By vector')
  if (length(scale)!=nVar) {
    scale=match.arg(scale,c('UV','Pareto','none'))
    cat('\nScale:',scale)
    if (scale=='UV') scale=TRUE else if (scale=='none') scale=FALSE else scale=apply(X,2,function(x) sqrt(sd(x)))
  } else cat('\nScale: By vector')
  if(!(is.logical(scale) | length(scale)==nVar)) stop('Error with scaling')
  X=scale(X,center,scale)
}
