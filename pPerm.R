#' Calculate permutation p-value of actual model performance vs null hypothesis distribution
#' 
#' `pPerm` will calculate the cumulative (1-tailed) probability of `actual` belonging to `h0`. 
#' Side is guessed by actual value compared to median(h0).
#' Test is performed on original data OR ranked for non-parametric statistics.
#' @param actual Actual model performance (e.g. misclassifications or Q2)
#' @param h0 Null hypothesis distribution from permutation test (same metric as `actual`)
#' @param side Smaller or greater than (automatically guessed if omitted) (Q2 is a "greater than" test, whereas misclassifications is "smaller than")
#' @param type Standard Student's t distribution ('t') or Student's t on rank-transformed data for nonparametric test ('non')
#' @return p-value
#' @export
pPerm=function(actual,h0,side=c('smaller','greater'),type=c('t','non')) {
  if(missing(type)) type='t'
  if(missing(side)) side=ifelse(actual<median(h0),'smaller','greater')
  if (type=='non') {
    rank=rank(c(actual,h0))
    actual=rank[1]
    h0=rank[-1]
  }
  p=pt((actual-mean(h0))/sd(h0),df=length(h0)-1)
  if (side=='greater') p=1-p
  return(p)
}
