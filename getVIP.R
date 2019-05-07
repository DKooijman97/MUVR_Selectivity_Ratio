#' Extract autoselected variables from MUVR model object
#'
#' @param MVObj 
#' @param model Which model to use ("min", "mid" (default), or "max")
#'
#' @return Data frame with order, name and average rank of variables (`order`, `name` & `rank`)
#' @export
getVIP=function(MVObj,model='mid') {
  nMod=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nVar=round(MVObj$nVar[nMod])
  VIPs=sort(MVObj$VIP[,nMod])[1:nVar]
  VIPs=data.frame(order=1:nVar,name=names(VIPs),rank=VIPs)
  VIPs$name=as.character(VIPs$name)
  return(VIPs)
}