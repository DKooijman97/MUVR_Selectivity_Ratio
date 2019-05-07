#' Plots variable importance ranking in MUVR object
#'
#' Regardless of MV core method, variables are sorted by rank, where lower is better. `plotVIP` produces boxplots of variable rankings for all model repetitions. 
#'
#' @param MVObj MUVR object
#' @param n Number of top ranking variables to plot (defaults to those selected by MUVR)
#' @param cut Optional value to cut length of variable names to `cut` number of characters
#' @param model Which model to choose ('min', 'mid' {default} or 'max')
#'
#' @return Barplot of variable rankings (lower is better)
#' @export
plotVIP=function(MVObj,n,model='mid',cut) {
  nModel=ifelse(model=='min',1,ifelse(model=='mid',2,3))
  nFeat=round(MVObj$nVar[nModel])
  if(missing(n)) n=nFeat
  VIP=MVObj$VIP[,nModel]
  VIPRep=MVObj$VIPPerRep[[nModel]]
  VIPRep=VIPRep[order(VIP),][1:n,]
  if(n>nFeat) {
    VIPRep=rbind(VIPRep[1:nFeat,],rep(NA,ncol(VIPRep)),VIPRep[(nFeat+1):n,])
    col=rep(c('yellow','grey'),c(nFeat,(n-nFeat+1)))
  } else col=NULL
  VIPRep=VIPRep[nrow(VIPRep):1,]
  col=rev(col)
  boxplot(t(VIPRep),horizontal=T,axes=F,col=col)
  axis(1)
  labels=rownames(VIPRep)
  if(!missing(cut)) labels=substring(labels,1,cut)
  axis(2,las=1,at=1:nrow(VIPRep),labels=labels)
  if (n>nFeat) abline(h=(n-nFeat+1))
  box(bty='o')
}
