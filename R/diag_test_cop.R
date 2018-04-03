#' Copula Diagnostic
#' 
#' \code{diagMAAR.cop} takes a preprocessed data matrix test the MAAR assumption
#' using a Gaussian copula approach - as outlined in our paper.
#' 
#' \code{diagMAAR.cop} is part of the diagnostic tools functions used for
#' diagnosing for the MAAR assumption. This particular one assumes that the
#' correlation between the (data, missingness indicators) can be modeled using a
#' Gaussian copula. This function utilizes the \code{sbgcop} package to sample
#' from the posterior distribution of the copula. The code then looks for 
#' conditional dependencies between the partially observed outcome and the 
#' different missingness indicators. 
#' 
#' @param prep A preprocessed S3 class that contains the data that is going to
#' be tested.
#' @param alpha A numeric value indicating the level of the test.
#' The default is set to 0.05.
#' @param nburn A numeric value indicating the number of burn in samples for the
#' MCMC (passed to sbgcop). The default is set to 500.
#' @param n.samples A numeric value indicating the number of samples to keep
#' from the MCMC (passed to sbgcop). The default is set to 2000.
#' @param verbose A logical variable indicating if sbgcop should print it's 
#' output. The default is set to FALSE.
#' @return diagMAAR A S3 object that contains: reject, a logical indicating if 
#' the test rejected; res, the results from the likelihood ratio test;
#'  which.reject, a vector indicating which variables were reject; method, a 
#' string indicating the diagnostic method used.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.mvn <- MANAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 1, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7,0.2, - 0.2), 
#'                                                 1, 5))
#' Y.cop <- prep.cop(obs.mvn)
#' diagMAAR.cop(Y.cop)
#' @family diagnostic
diagMAAR.cop <- function(
  prep, alpha = 0.05, nburn = 500, n.samples = 2000, verbose = FALSE) {
  nvmis <- prep$nvar.miss
  # Obtain samples from the posterior 
  OUT <- sbgcop::sbgcop.mcmc(prep$dt, nsamp = n.samples, verb = verbose)
  # Extract the correlation 
  Tstat <- array(rep(0, ((nvmis ^ 2) ^ 2 * (dim(OUT$C.psamp)[3] - nburn))),
                 dim = c(nvmis * 2, nvmis * 2, (dim(OUT$C.psamp)[3] - nburn)))
  # Compute the conditional posterior distribution.
  for (k in nburn:(dim(OUT$C.psamp)[3])) {
    Tstat[, , k - nburn] <- cond_cov(OUT$C.psamp[, , k], nvmis * 2)
  }
  # Store the results in a array 
  results <- array(0, dim = c(nvmis, nvmis, 3))
  # If there are more than one tests use a Bonferroni correction
  if (nvmis > 1) {
    BFcorrect <- nvmis ^ 2 - nvmis
  }else {
    BFcorrect <- 1
  }
  for (k in 1:nvmis) {
    for (l in 1:nvmis) {
      results[k,l,] <- extract_stat(
          Tstat[k, l + nvmis, ], (alpha / BFcorrect))
    }
  }
  # Check if any of the intervals fall outside of zero 
  mat1 <- results[, , 1]
  mat2 <- results[, , 3]
  # If there are more than 1 test we remove the diagonals as they can not be
  # estimated. 
  if (nvmis > 1) {
    diag(mat1) <- 0
    diag(mat2) <- 0
  }
  which.reject <- apply(rbind(apply(mat1 < 0 & mat2 < 0, 2, any), 
                      apply(mat1 > 0 & mat2 > 0, 2, any)), 2, any)
  reject <-  any(mat1 < 0 & mat2 < 0 ) || any(mat1 > 0 & mat2 > 0)
  diagMAAR <- list(reject = reject, results = results, 
                   which.reject = which.reject, method = "cop")
  class(diagMAAR) <- "diagMAAR"
  return(diagMAAR)
}

#' extract_stat
#' 
#' This function is used to extract the mean, medium, upper and lower quantiles
#' (as specified by alpha) of the input vector.
#' 
#' @param Vec A numeric matrix, usually a covariance matrix
#' @param alpha An numeric value specifying the upper and lower quantiles
#' @param med A logical value, if TRUE the function returns the medium instead 
#' of the  the mean.
extract_stat <- function(Vec, alpha = 0.05, med = FALSE) {
  alpha <- alpha / 2
  if (med) {
    res <- rep(0,4)
    res[1] <- mean(Vec)
    res[2:4] <- quantile(Vec, c(alpha, .5, (1 - alpha)))
  }else {
    res <- quantile(Vec, c(alpha, .5, (1 - alpha)))
  }
  return(res)
}

#' cond_cov
#' 
#' This function is used to compute the conditional covariance when using the
#' \code{diagMAAR.cop} method.
#' 
#' @param Met A numeric matrix, usually a covariance matrix
#' @param nv An integer indicating the location of what is to be inverted.
cond_cov<-function(Met, nv) {
  # Compute the conditional distribution
  n <- dim(Met)[1]
  Vcond <- Met[1:nv, 1:nv] - Met[1:nv, (nv + 1):n] %*%
    solve(Met[(nv + 1):n, (nv + 1):n]) %*% Met[(nv + 1):n, 1:(nv)]
  return(Vcond)
}




