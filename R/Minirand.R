#' Minimization randomization to k treatment groups
#' 
#' The function is used to generate treatment assignment by minimization algorithms.
#' @param covmat matrix or data frame of covariate factors
#' @param j the jth subject in the randomization sequence
#' @param covwt vector of weights of the covaraite factors
#' @param ratio vector of randomization ratios for each treatment
#' @param ntrt numeric number of treatment groups
#' @param trtseq vector of a sequence of treatment groups
#' @param method the method or algorithm for the minimization randomization
#' @param result the treatment assignments in subjetcs achieved so far
#' @param p the high probability for new assignment
#' @return treatment assignment for the jth subject
#' @references Pocock and Simon (1975), Sequential Treatment Assignment with Balancing for Prognostic
#'  Factors in the Controlled Clinical Trial. Biometrics; 103-115.
#' @references Jin, Polis, and Hartzel (2019). "Algorithms for minimization 
#' randomization and the implementation with an R package". Communications in Statistics-Simulation 
#' and Computation; May 2019. 
#' @examples ntrt <- 3
#' nsample <- 120
#' trtseq <- c(1, 2, 3)
#' ratio <- c(2, 2, 1)
#' c1 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.4, 0.6)) 
#' c2 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.3, 0.7))
#' c3 <- sample(c(2, 1, 0), nsample, replace = TRUE, prob = c(0.33, 0.2, 0.5)) 
#' c4 <- sample(seq(1, 0), nsample, replace = TRUE, prob = c(0.33, 0.67)) 
#' covmat <- cbind(c1, c2, c3, c4) # generate the matrix of covariate factors for the subjects
#' # label of the covariates 
#' colnames(covmat) = c("Gender", "Age", "Hypertension", "Use of Antibiotics") 
#' covwt <- c(1/4, 1/4, 1/4, 1/4) #equal weights
#' res <- rep(100, nsample) # result is the treatment needed from minimization method
#' #gernerate treatment assignment for the 1st subject
#' res[1] = sample(trtseq, 1, replace = TRUE, prob = ratio/sum(ratio)) 
#' for (j in 2:nsample)
#' {
#' # get treatment assignment sequentiall for all subjects
#' res[j] <- Minirand(covmat=covmat, j, covwt=covwt, ratio=ratio, 
#' ntrt=ntrt, trtseq=trtseq, method="Range", result=res, p = 0.9)
#' }
#' trt1 <- res
#' #Display the number of randomized subjects at covariate factors
#' balance1 <- randbalance(trt1, covmat, ntrt, trtseq) 
#' balance1
#' totimbal(trt = trt1, covmat = covmat, covwt = covwt, 
#' ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range")
#' @export

Minirand <- function(covmat = covmat, j, covwt = covwt, ratio = ratio, ntrt = ntrt, trtseq = trtseq, method = "Range", result = res, p)
{
  if (j > 1) {
    matchx = apply(covmat[1:(j - 1), , drop = FALSE], 1, 
                   function(x, xrow) {
                     as.numeric(x == xrow)
                   }, covmat[j, ])    
    n_matchx <- matrix(0, ncol(covmat), ntrt)
    for (k in 1:ntrt)
    {
      n_matchx[,k] <- apply(as.matrix(matchx[, result[1:(j-1)]==trtseq[k]]), 1, sum)
    }   
    
    if (method == "Range")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        range_level <- apply(num_level, 1, range)
        imb_margin <- range_level[2, ] - range_level[1, ]
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }  
    
    if (method == "SD")
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[, i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        sd_level <- apply(num_level, 1, sd)
        imb_margin <- sd_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    if (method == "Var") 
    {
      imbalance <- rep(0, ntrt)
      for (i in 1:ntrt)
      {
        temp <- n_matchx
        temp[,i] <- temp[, i]+1 
        num_level <- temp %*% diag(1/ratio)   
        var_level <- apply(num_level, 1, var)
        imb_margin <- var_level
        imbalance[i] <- sum(covwt %*% imb_margin)
      }
    }
    
    trt.mini <- trtseq[imbalance == min(imbalance)]  
    trt.highprob <- trt.mini 
    trt.lowprob <- trtseq[-trt.mini]
    res <- ifelse(length(trt.highprob) < ntrt, sample(c(trt.highprob, trt.lowprob), 1, replace = TRUE, prob = c(rep(p/length(trt.highprob), length(trt.highprob)),
                                                                                                                rep((1-p)/length(trt.lowprob), length(trt.lowprob)))), sample(trtseq, 1, replace = TRUE, prob = rep(1/ntrt, ntrt)))
  }
  return(res)
}