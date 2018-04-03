
#' MAAR diagnostics
#' 
#' \code{diagMAAR} takes a preprocessed data matrix test the MAAR assumption.
#' 
#' \code{diagMAAR} implements three different diagnostic test for the missing
#' always at random assumption: Gaussian copula test (cop), directly testing
#' a postulated missingness mechanism (dtmm), and comparison of conditional 
#' means (ccm). Based on simulation studies the dtmm method should be avoided
#' for small sample sized (< 100). For more details about each of the different
#' diagnostic test see diagMAAR.test_name.
#' 
#' @param dt A data.frame with missing values, a matrix with missing values,
#' a mids class (the output from mice), or a prep class (the output from) prep.
#' @param daigtest The test to be implemented - pick from Gaussian copula test 
#' (cop), directly testing a postulated missingness mechanism (dtmm), and 
#' comparison of conditional means. Leave blank if dt is of class prep.
#' @param alpha A numeric value indicating the level of the test.
#' The default is set to 0.05.
#' @param ... Other arguments to be passed on to the specific diagnostic method
#' @param seed A numeric value used to set a seed, leave blank for NULL.
#' @return A S3 object that contains: reject, a logical indicating if 
#' the test rejected; res, the results from the likelihood ratio test;
#'  which.reject, a vector indicating which variables were reject; method, a 
#' string indicating the diagnostic method used.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.nvm <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 4, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
#' diagMAAR.dtmm(obs.mvn, "ccm")
#' @family diagnostic
#' @importFrom sbgcop sbgcop.mcmc
diagMAAR <- function(dt, daigtest = NULL, alpha = 0.05, ..., seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  if (class(dt) == "mids") {
    dt <- prep(dt$data, daigtest = daigtest)
  } else if ("data.frame" %in% class(dt) | "matrix" %in% class(dt)) {
    dt <- prep(data, daigtest = daigtest)
  }
  if (class(dt) == "prep") {
    if (dt$daigtest == "cop") {
      diagMAAR <- diagMAAR.cop(dt, ...)
    } else if (dt$daigtest == "dtmm") {
      diagMAAR <- diagMAAR.dtmm(dt, ...)
    }  else if (dt$daigtest == "ccm") {
      diagMAAR <- diagMAAR.ccm(dt, ...)
    } else {
      cat("Invalid diagnostic test specified.")
      cat("Please use cop, dtmm, or ccm only.")
      diagMAAR <- NULL
    }
  }
  return(diagMAAR)
}

#' @export
print.diagMAAR <- function(diagMAAR, ...) {
  cat("Diagnostic Test: ")
  print(diagMAAR$daigtest)
  cat("\nReject: ")
  print(diagMAAR$reject)
  cat("\nReject by variable:")
  print(diagMAAR$which.reject)
}








