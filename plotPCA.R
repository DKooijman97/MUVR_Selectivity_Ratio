#' PCA score plot from prcomp()
#'
#' Customised PCA score plots with the possibility to choose PCs, exporting to png and the possibility to add color or different plotting symbols according to variable.
#' @param pca A `prcomp` object
#' @param PC1 Principal component on x-axis
#' @param PC2 Principal component on y-axis
#' @param file If specified provides the name of a png export file. Otherwise normal plot.
#' @param colVar Continuous variable for coloring observations (40 cuts)
#' @param symbVar Categorical/discrete variable for multiple plot symbols
#' @param main If provided provides a main title of the plot
#'
#' @return A PCA score plot. Exported as png if `file` specified in function call.
#' @export
plotPCA=function(pca,PC1=1,PC2=2,file,colVar,symbVar,main='') {
  if (missing(colVar)) {
    col=1
    colLeg=FALSE
  } else {
    cols=colorRampPalette(c('blue','yellow','red'))(40)
    col=cols[cut(colVar,40)]
    colLeg=TRUE
  }
  if (missing(symbVar)) {
    pch=1
    symbLeg=FALSE
  } else {
    symbVar=factor(symbVar)
    symbs=c(1,2,0,6)
    nSymb=length(levels(symbVar))
    if (nSymb>4) {
      symbs=c(symbs,9:(nSymb+4))
    }
    pch=symbs[symbVar]
    symbLeg=TRUE
  }
  plotPNG=ifelse(missing(file),FALSE,TRUE)
  if (plotPNG) png(filename=file,width=1024,height=1024,pointsize=36)
  if (main=='') par(mar=c(4,4,0,0)+.5) else par(mar=c(4,4,2,0)+.5)
  pcVar=summary(pca)$importance[2,]
  xlab=paste('PC',PC1,' (R2X=',signif(pcVar[PC1],3),')',sep='')
  ylab=paste('PC',PC2,' (R2X=',signif(pcVar[PC2],3),')',sep='')
  plot(pca$x[,c(PC1,PC2)],main=main,col=col,xlab=xlab,ylab=ylab,pch=pch) # scoreplot
  if (colLeg) legend('topleft',col=c('blue','yellow','red'),legend=c('low','mid','high'),pch=1)
  if (symbLeg) legend('bottomleft',col=1,legend=levels(symbVar),pch=symbs[1:length(levels(symbVar))])
  abline(h=0,lty=2)
  abline(v=0,lty=2)
  if (plotPNG) dev.off()
}
