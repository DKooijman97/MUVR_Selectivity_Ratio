#' Extract PLS Sratio values
#' Authors: Dennis Kooijman
#' Date: 
#' 
#' Calculates of selectivity ratio to rank variables, as alternitive
#' to the VIP approach. Based on the orginal vip module of the MUVR algorithm (as downloaded on 6-5-2019). 
#'
#' Saves the Sratio in a vip object because this object plays crucial role 
#' in the other parts of the program
#' 
#' @param object pls(da) object
#'
#' @return vip object (in which the Sratio is stored)
#' @export

Sratio <- function(object) {
  #-- initialisation of needed matrices --#
  #b = object$mat.c                   # Matrix with coefficients (nrow = number of variables, ncol = number of PCs)
  H = object$ncomp                    # H is number of latent variables (componements)      
  X = object$X                        # The preprocessed X-data
  W = 
  scores  = object$variates$X         # Scores (t)
  loadings = object$loadings$X        # Loadings (p)
  
  #Getting number of samples and variables
  n = nrow(X)  
  nVariables = ncol(X)
  
  # Calculation of the normalized PLS regression vector  
  #  for (i in 1:h) { 
  #      b[,h] = b[,h]/sqrt(sum(b[,h]^2))
  #  }
  
  # Initialization of Sratio matrix  
  Sratio = matrix(0, nrow = nVariables, ncol = 1 )
  
  # Calculation of Sratio per variable per sample for every latent variable
  SS_explained = (scores%*%t(loadings))^2
  
  #SS_explained = colSums(SS_explained*SS_explained)
  
 
  SS_unexplained = colSums(X-SS_explained); 
   
  
  for (i in 1:nVariables) {
     Sratio[i]=SS_explained[i]/SS_unexplained[i]
     
  }
  
  
  # Returning Sratio in VIP object 
  rownames(Sratio) = rownames(loadings)
  colnames(Sratio)= "Comp 1"
  return(invisible(Sratio))
}
