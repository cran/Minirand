#' Displays the number of randomized subjects at each level for all covariate factors.
#' 
#' The fuction to cound the number of randomized subjects at each level for all covariate factors
#' @param trt treatment sequence for all the randomized subjects
#' @param covmat matrix or data frame of covariate factors
#' @param ntrt numeric number of treatment groups
#' @param trtseq vector of a sequence of treatment groups
#' @return the number of randomized subjects at each level for all covariate factors
#' @export

randbalance <- function(trt, covmat, ntrt, trtseq)
{
  balance <- vector(length = ncol(covmat), "list")
  names(balance) = colnames(covmat)
  for (i in 1:ncol(covmat))
  {
    
    balance[[i]] <- table(trt, covmat[, i])
    
  }
  return(balance)
}