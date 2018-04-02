#' Sample a Multivariate Normal
#' 
#' This function creates N samples from a multivariate normal distribution
#' with mean equal to zero and pairwise-correlation rho.
#' 
#' This is a basic functions that samples a standard normal and uses the
#' Choleski decomposition to create correlated samples.
#' 
#' @param nvar An integer specifying the dimension of the multivariate normal
#'         distribution.
#' @param rho The correlation matrix, if a single number is given then the
#'        correlation is assumed to be constant between all of the variable.
#' @param N An integer specifying the number of samples.
#' @param seed The seed used to generate that data set.seed.
#' @keywords Multivariate Normal Sample
#' @return samples.mvn A matrix of dimensions equal to (nvar, N) containing the
#'         samples from a MVN.
#' @export
#' @examples sample_mvn(5, 0.3, 100)
#' @family data generating functions
sample_mvn <- function(nvar, rho, N, seed = 0) {
  set.seed(seed)
  if (is.null(dim(rho))) {
    # If rho is a number compute the correlation matrix.
    rho <- matrix(rho, nvar, nvar)
    diag(rho) <- 1
  }
  # Compute the Choleski decomposition of rho. 
  rho.chol <- chol(rho)
  # Generate the multivariate data set
  samples.mvn <- matrix(rnorm(N * nvar), ncol = nvar) %*% rho.chol
  return(samples.mvn)
}
