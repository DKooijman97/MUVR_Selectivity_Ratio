#' PLS-DA
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#'
#' @param x 
#' @param y 
#' @param ncomp 
#' @param max.iter 
#' @param tol 
#' @param near.zero.var 
#' @param scale 
#'
#' @return plsda object       # used in vip module, DennisK
#' @export
plsda <- function(x, y, ncomp = 2, max.iter = 500, tol = 1e-06, near.zero.var = TRUE, scale = TRUE) {
  y = as.factor(y)	

  n <- length(y)
  groups <- sort(unique(y))
  levels =  levels(y)### Add levels
  cgroups <- as.character(groups)
  groups <- as.numeric(factor(cgroups, levels = unique(cgroups)))
  classification <- as.numeric(factor(as.character(y), levels = unique(cgroups)))
  k <- length(groups)
  nam <- levels(groups)
  ind.mat <- matrix(0, n, k, dimnames = c(names(classification), nam))
  for (j in 1:k) ind.mat[classification == groups[j], j] <- 1
  attr(ind.mat, "levels") = levels
  
  result = MUVR::pls(x, ind.mat, ncomp = ncomp, max.iter = max.iter, tol = tol,near.zero.var = near.zero.var, scale = scale)
  
  result$ind.mat = ind.mat
  result$names$Y = levels(y)
  class(result) = "plsdaMUVR"
  return(invisible(result))	
}
