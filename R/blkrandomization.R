#' Blocked randomization
#' 
#' The fuction is used to generate treatment assignments based on blocked randomization. 
#' @param n numeric number of subjects who will be randomized
#' @param blocksize numeric value of block size used for blocked randomization
#' @param block vector of treatment blocks used for blocked randomization
#' @return trt a sequence of treatment assignments
#' @examples blocksize <- 4
#' block <- c(1, 2, 3, 4) # treatment 1, 2, 3, 4
#' n <- 35
#' blkrandomization(n, blocksize, block)
#' @export

 
blkrandomization <- function(n, blocksize, block)
{
  trt <- rep(100, n)
  number_of_blocks <- as.integer(n/blocksize)
  if (number_of_blocks > 0)
  {
    random <- NULL
    for (i in 1:number_of_blocks)
    {
      random <- c(random, sample(block, replace = FALSE))
    }
    trt[1:(number_of_blocks * blocksize)] <- random
    if (n > (number_of_blocks * blocksize))
    {
      trt[(number_of_blocks * blocksize+1):n] <- sample(block, length((number_of_blocks * blocksize + 1):n), replace = FALSE)
    }
  }
  trt
}
