#' Directly testing a postulated missingness mechanism.
#' 
#' \code{diagMAAR.dtmm} takes a preprocessed data matrix test the MAAR 
#' assumption using a logistic regression as the response propensity model.
#' 
#' \code{diagMAAR.dtmm} is part of the diagnostic tools functions used for
#' diagnosing for the MAAR assumption. This function tests if the response 
#' propensity in one variable depends on partially observed outcome 
#' variables. To perform the likelihood ratio test, the function first imputes
#' the missing values - using default mice settings.
#' 
#' Note: In simulation studies this test has low power for sample sizes < 100.
#' 
#' @param prep A preprocessed S3 class that contains the data that is going to
#' be tested.
#' @param alpha A numeric value indicating the level of the test.
#' The default is set to 0.05.
#' @keywords MAAAR testing
#' @return diagMAAR A S3 object that contains: reject, a logical indicating if 
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
#' Y.dtmm <- prep.dtmm(obs.mvn)
#' diagMAAR.dtmm(Y.dtmm)
#' @family diagnostic
#' @import mice

diagMAAR.dtmm <- function(prep, alpha = 0.05, verbose = FALSE) {
  # Run a quick prediction to set up the imputation matrix
  pred <- mice::quickpred(prep$dt)
  # Make sure we don't use R to predict 
  pred[, (prep$nvar + 1):(prep$nvar + prep$nvar.miss)] <- 0
  # Run the imputations using a Gaussian model. 
  imp <- mice::mice(prep$dt, pred = pred, m = 10, maxit = 10, print = FALSE)
  res <- NULL
  for (kk in 1:prep$nvar.miss) {
    if (verbose) {
      cat(kk)
      cat("\r")
    }
    fit_R <- with(imp, glm(as.formula(paste0("R", as.character(kk), " ~ ", 
                          paste(sapply(setdiff(c(1:prep$nvar), prep$locmis), 
                                 function(ll) paste0("Y", as.character(ll))),
                                collapse = " + "), " + ", 
                          paste(sapply( prep$locmis, 
                                 function(ll) paste0("Y", as.character(ll))),
                                collapse = " + "))), family = binomial()))
    fit_Rs <- with(imp, glm(as.formula(paste0("R", as.character(kk), " ~ ", 
                          paste(sapply(setdiff(c(1:prep$nvar), prep$locmis), 
                                 function(ll) paste0("Y", as.character(ll))),
                                collapse = " + "))), family = binomial()))
    res <- c(res, 
      mice::pool.compare(fit_R, fit_Rs,method = "likelihood",data = imp)$pvalue)
  }
  which.reject <- res < alpha / prep$nvar.miss
  reject <- any(which.reject)
  diagMAAR <- list(reject = reject, results = res, 
                   which.reject = which.reject, method = "dtmm")
  class(diagMAAR) <- "diagMAAR"
  return(diagMAAR)
}
 