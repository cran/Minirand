#' Calculates the total imbalance measured by minimization algorithms.
#' 
#' The function to calculates the total imbalance measured by minimization algorithms
#' @param trt treatment sequence for all the randomized subjects
#' @param covmat matrix or data frame of covariate factors
#' @param covwt vector of weights of the covaraite factors
#' @param ratio vector of randomization ratios for each treatment
#' @param ntrt numeric number of treatment groups
#' @param trtseq vector of a sequence of treatment groups
#' @param method the method or algorithm for the minimization randomization
#' @return total imbalance
#' @export


totimbal <- function (trt = trt, covmat = covmat, covwt = covwt, ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
{
  balance <- randbalance(trt, covmat, ntrt, trtseq)
  if (method == "Range")
  {
    imbalance<-rep(0, ncol(covmat)) 
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp) %*% diag(1/ratio)   
      range_level <- apply(num_level, 1, range)
      imb_margin <- range_level[2,] - range_level[1,]
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }  
  
  if (method == "SD")
  {
    imbalance <- rep(0, ncol(covmat))
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp)%*%diag(1/ratio)   
      sd_level <- apply(num_level, 1, sd)
      imb_margin <- sd_level
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }
  
  if (method == "Var")
  {
    imbalance <- rep(0, ncol(covmat))
    for (i in 1:ncol(covmat))
    {
      temp <- balance[[i]]
      num_level <- t(temp) %*% diag(1/ratio)   
      var_level <- apply(num_level, 1, var)
      imb_margin <- var_level
      imbalance[i] <- sum(imb_margin)
    }
    total_imb <- sum(imbalance*covwt)
  }
  return(total_imb)
}