#' PLS biplot
#'
#' Makes a biplot of a fitted object (e.g. from a MUVR with PLS core)
#'
#' @param fit A PLS fit (e.g. from MVObject$Fit[[2]])
#' @param comps Which components to plot
#' @param xCol (Optionsal) Continuous vector for grey scale gradient of observation (sample) color (e.g. Y vector in regression analysis)
#' @param labPlSc Boolean to plot observation (sample) names (defaults to TRUE)
#' @param labs (Optional) Labels names
#' @param vars Which variables to plot (names in rownames(loadings))
#' @param labPlLo Boolean to plot variable names (defaults to TRUE)
#' @param pchSc Plotting character for observation scores
#' @param colSc Colors for observation scores (only if xCol omitted)
#' @param colLo Colors for variable loadings (defaults to red)
#' @param supLeg Boolean for whether to suppress legends
#'
#' @return A PLS biplot
#' @export
biplotPLS=function(fit,comps=1:2,xCol,labPlSc=TRUE,labs,vars,labPlLo=TRUE,pchSc=16,colSc,colLo=2,supLeg=FALSE) {
  cex=par()$cex
  par(mar=c(4,4,4,4)+.5)
  scores=fit$variates$X[,comps]
  loads=fit$loadings$X[,comps]
  if(missing(vars)) vars=rownames(loads)
  loads=loads[rownames(loads)%in%vars,]
  nSamp=nrow(scores)
  nVar=nrow(loads)
  if(missing(xCol)) {
    if (missing(colSc)) {
      colSc=rep(1,nSamp)
      legPlot=FALSE
    } else {
      colScLeg=colSc
      colSc=as.factor(colSc)
      legPlot=TRUE
    }
  } else {
    x.col=10+round(85*((max(xCol)-xCol)/(diff(range(xCol)))))
    colSc=paste("gray",x.col,sep="")
    legPlot=TRUE
  }
  if(supLeg) legPlot=FALSE
  rLo=max(abs(loads))
  rLo=1.1*c(-rLo,rLo)
  plot(loads,xlim=rLo,ylim=rLo,type='n',xlab='',ylab='',main='',axes=F)
  box(bty='o')
  axis(3)
  axis(4,las=1)
  mtext('Loadings',3,line=3,cex=cex)
  mtext('Loadings',4,line=3,cex=cex)
  abline(h=0,v=0,lty=2,col='grey')
  arrows(rep(0,nrow(loads)),rep(0,nrow(loads)),loads[,1],loads[,2],col=colLo)
  if(labPlLo) text(1.1*loads,rownames(loads),cex=.7*cex,font=3)
  par(new=T)
  rSc=max(abs(scores))
  rSc=c(-rSc,rSc)
  plot(scores,col=colSc,pch=pchSc,ylim=rSc,xlim=rSc,bty='l',xlab=paste('Component',comps[1],'Scores'),ylab=paste('Component',comps[2],'Scores'),axes=F)
  axis(1)
  axis(2,las=1)
  if (labPlSc) {
    if (missing(labs)) labs=rownames(scores)
    # cat(labs)
    text(scores,as.character(labs),pos=3)
  }
  if (legPlot) {
    if (missing(xCol)) {
      legend('topright',legend=unique(colScLeg),pch=unique(pchSc),col=unique(colSc),bty='n')
    } else {
      whUnik=!duplicated(xCol)
      unik=xCol[whUnik]
      cols=colSc[whUnik][order(unik)]
      unik=sort(unik)
      if (length(unik)>10) {
        k=(length(unik)-1)/5
        n=1+k*0:5
        cols=cols[n]
        unik=unik[n]
      }
      legend('topright',legend=signif(unik,3),fill=cols,bty='n')
    }
  }
}
