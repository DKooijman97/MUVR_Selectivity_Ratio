#' PLS regression 
#' 
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#'
#' @param x             #X data, (input in plsinner == Xtrain) DennisK 
#' @param y             #Y data (response)                     DennisK 
#' @param ncomp         #Number of latent variables            DennisK 
#' @param max.iter      #Max. number of iteration              DennisK
#' @param tol           #Tolerance                             DennisK
#' @param near.zero.var 
#' @param scale 
#'
#' @return pls object
#' @export
pls <- function(x, y, ncomp = 2, max.iter = 500, tol = 1e-06, near.zero.var = TRUE, scale = TRUE) {
  y = as.matrix(y)
  x = as.matrix(x)
  # Remove variables with near zero variance 
  if(near.zero.var) {
    nzv = MUVR::nearZeroVar(x)
    if (length(nzv$Position > 0)) {
      warning("Zero- or near-zero variance predictors.\nReset predictors matrix to not near-zero variance predictors.\nSee $nzv for problematic predictors.")
      x = x[, -nzv$Position,drop=FALSE]
      if(ncol(x)==0) stop("No more predictors after Near Zero Var has been applied!")
    }
  }
  n = nrow(x)
  p = ncol(x)
  q = ncol(y)
  
  # Names
  x.names = colnames(x)
  if (is.null(x.names)) {
    x.names = paste("X", 1:p, sep = "")
    colnames(x)=x.names
  }
  if (dim(y)[2] == 1) {
    y.names = "Y"
  } else {
    y.names = colnames(y)
    if (is.null(y.names)) {
      y.names = paste("Y", 1:q, sep = "")
      colnames(y)=y.names
    }
  }
  ind.names = rownames(x)
  if (is.null(ind.names)) {
    ind.names = rownames(y)
    rownames(x) = ind.names
  }
  if (is.null(ind.names)) {
    ind.names = 1:n
    rownames(x) = rownames(y) = ind.names
  }	
  
  # Center and scale indata
  if (scale) x.temp = x = scale(x, center = TRUE, scale = TRUE) else x.temp=x
  y.temp = y = scale(y, center = TRUE, scale = TRUE) 
  
  # Allocate matrices
  mat.t = mat.u = matrix(nrow = n, ncol = ncomp)
  mat.a = mat.c = matrix(nrow = p, ncol = ncomp)
  mat.b = mat.d = matrix(nrow = q, ncol = ncomp)
  
  # Iterate pls components h
  iter=NULL
  for (h in 1:ncomp) {
    #-- initialisation --#
    M = crossprod(x.temp, y.temp)
    svd.M = svd(M, nu = 1, nv = 1)
    a.old = svd.M$u
    b.old = svd.M$v
    #-- latent variables --#
    t = x.temp %*% a.old / drop(crossprod(a.old))
    u = y.temp %*% b.old / drop(crossprod(b.old))
    iterh = 1
    #-- convergence of a  --#
    repeat{
      a = t(x.temp) %*% u 
      a = a / drop(sqrt(crossprod(a)))
      t = x.temp %*% a / drop(crossprod(a))
      b = t(y.temp) %*% t 
      b = b / drop(sqrt(crossprod(b)))
      u = y.temp %*% b / drop(crossprod(b))
      if (crossprod(a - a.old) < tol) break
      if (iterh == max.iter) break
      a.old = a
      b.old = b
      iterh = iterh + 1
    }
    #-- deflation --#
    c = crossprod(x.temp, t) / drop(crossprod(t))
    x.temp = x.temp - t %*% t(c)   
    #-- mode regression --#
    d = crossprod(y.temp, t) / drop(crossprod(t))
    y.temp = y.temp - t %*% t(d)
    
    mat.t[, h] = t
    mat.u[, h] = u
    mat.a[, h] = a
    mat.b[, h] = b
    mat.c[, h] = c
    mat.d[, h] = d
    iter=c(iter,iterh) #save the number of iteration per component
  } 
  #-- valeurs sortantes --#     DennisK, the creation of the plsda object used in VIP
  rownames(mat.a) = rownames(mat.c) = x.names
  rownames(mat.b) = y.names
  rownames(mat.t) = rownames(mat.u) = ind.names
  comp = paste("comp", 1:ncomp)
  colnames(mat.t) = colnames(mat.u) = comp
  colnames(mat.a) = colnames(mat.b) = colnames(mat.c) = comp 
  result = list(X = x, Y = y, ncomp = ncomp, mat.c = mat.c, variates = list(X = mat.t, Y = mat.u), 
                loadings = list(X = mat.a, Y = mat.b), tol = tol, max.iter = max.iter, iter=iter)
  if (near.zero.var == TRUE) result$nzv = nzv
  class(result) = "plsMUVR"
  return(invisible(result))
}

