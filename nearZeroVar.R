#' Identify variables with nea zero variance
#'
#' Adapted and stripped down from mixOmics v 5.2.0 (https://cran.r-project.org/web/packages/mixOmics/)
#' @param x 
#' @param freqCut 
#' @param uniqueCut 
#'
#' @return nzv object
#' @export
nearZeroVar <- function (x, freqCut = 95/5, uniqueCut = 10) {
  if (is.vector(x)) x = matrix(x, ncol = 1)
  
  freqRatio = apply(x, 2, function(data) {
    if (length(unique(data)) == length(data)){ # No duplicate
      return(1)
    } else if (length(unique(data)) == 1) { # Same value
      return(0)
    } else {
      t = table(data)
      return(max(t, na.rm = TRUE)/max(t[-which.max(t)], na.rm = TRUE))
    }
  })
  
  lunique = apply(x, 2, function(data) length(unique(data)))
  percentUnique = 100 * lunique/nrow(x)
  zeroVar = (lunique == 1) | apply(x, 2, function(data) all(is.na(data)))
  
  out = list()
  out$Position = which((freqRatio > freqCut & percentUnique <= uniqueCut) | zeroVar)
  names(out$Position) = NULL
  out$Metrics = data.frame(freqRatio = freqRatio, percentUnique = percentUnique)
  out$Metrics = out$Metrics[out$Position, ]
  out
}
