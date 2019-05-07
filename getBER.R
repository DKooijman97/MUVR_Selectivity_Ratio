#' get Balanced Error Rate from classification analysis
#'
#' @param actual Vector of actual classifications of samples
#' @param predicted Vector of predicted classifications of samples
#'
#' @return Balanced Error Rate (BER)
#' @export
#'
getBER=function (actual, predicted) 
{
  if (length(actual) != length(predicted)) stop ("Mismatch in length of arguments")
  if (!is.factor(actual)) actual = factor(actual)
  levs = levels(actual)
  nlevs = length(levs)
  confMat = matrix(0, nrow=nlevs, ncol=nlevs + 1)
  rownames(confMat) = levs
  colnames(confMat) = paste0("pred.", c(levs, "NA"))
  for (i in 1:nlevs) {
    whLev.i = which(actual == levs[i])
    for (j in 1:nlevs) confMat[i, j] = sum(predicted[whLev.i] == levs[j], na.rm = TRUE)
    confMat[i, nlevs + 1] = sum(is.na(predicted[whLev.i]))
  }
  if (sum(is.na(predicted)) == 0) confMat = confMat[, -(nlevs + 1)]
  confMat.wrong = confMat
  diag(confMat.wrong) = 0
  BER = sum(apply(confMat.wrong, 1, sum, na.rm = TRUE)/apply(confMat, 1, sum, na.rm = TRUE), na.rm = TRUE)/nlevs
  return(BER)
}
