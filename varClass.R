#' Report variables belonging to different classes
#' 
#' Reports names and numbers of variables: all as well as Optimal (min model), redundant (from min up to max) and noisy (the rest)
#' @param MVObject A MUVR object
#'
#' @return A list with names and numbers of variables: all as well as Optimal (Corresponding to 'min' or minial-optimal model), Redundant (from 'min' up to 'max' or all-relevant ) and Noisy (the rest)
#' @export
#'
varClass=function(MVObject) {
  nVarO=round(MVObject$nVar[1])
  nVarOR=round(MVObject$nVar[3])
  O=names(sort(MVObject$VIP[,1])[1:nVarO])
  OR=names(sort(MVObject$VIP[,3])[1:nVarOR])
  R=OR[!OR%in%O]
  ALL=rownames(MVObject$VIP)
  N=ALL[!ALL%in%c(O,R)]
  numbers=c(length(ALL),length(O),length(R),length(N))
  names(numbers)=c('All','Optimal','Redundant','Noisy')
  return(list(ALL=ALL,Optimal=O,Redundant=R,Noisy=N,numbers=numbers))
}