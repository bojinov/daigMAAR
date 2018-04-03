#' Comparison of Conditional Means
#' 
#' \code{diagMAAR.ccm} takes a preprocessed data matrix test the MAAR assumption
#' by comparing the conditional means.
#' 
#' \code{diagMAAR.ccm} is part of the diagnostic tools functions used for
#' diagnosing for the MAAR assumption. This function looks for a difference in
#' the conditional means for a variable. This test assume that the data are 
#' MAAR, the columns of the missingness indicators are mutually
#' conditionally independent given the outcome matrix, and the rows are
#' exchangeable.  Consider three variables; \eqn{Y_{1}} and \eqn{Y_{2}} can have
#' missing values, and \eqn{Y_{3}} is always fully observed. Then conditional 
#' expectation of \eqn{Y_{1}} given \eqn{Y_{3}} is the same for the two 
#' partitions induced by \eqn{R_{2}},
#' \deqn{E[Y_{i,1} | Y_{i,3} = y_{i,3}, R_{i,2} = 0] 
#'       = E[Y_{i,1} | Y_{i,3} = y_{i,3}, R_{i,2} = 1] 
#'       = E[Y_{i,1} | Y_{i,3} = y_{i,3}].}
#' Then, we can test if \eqn{Y_1} depends on \eqn{R_2} given \eqn{Y_3}. 
#' Similarly we can test if \eqn{Y_2} depends on \eqn{R_1} given \eqn{Y_3}. 
#' The test is assumes that the variables with missing values are Gaussian
#' and performs a likelihood ratio test between the model that includes only the
#' fully observed variables and a model that has an interaction between the
#' fully observed variables and the missingness indicators. 
#' Future version will include other models. 
#' 
#' 
#' @param prep A preprocessed S3 class that contains the data that is going to
#' be tested.
#' @param alpha A numeric value indicating the level of the test.
#' The default is set to 0.05.
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
#' Y.ccm <- prep.ccm(obs.mvn)
#' diagMAAR.ccm(Y.ccm)
#' @family diagnostic
diagMAAR.ccm <- function(prep, alpha = 0.05){
  pval <- diag(prep$nvar.miss)
  for (jj in 1:prep$nvar.miss) {
    # Create a new data set that doesn't include the missing variables or R
    dts.temp <- data.frame(out = prep$dt[, prep$locmis[jj]], 
                                 prep$dt[, - prep$locmis])
    # Fit the smaller model
    fit.small <- lm(out ~ ., dts.temp)
    for(kk in 1:prep$nvar.miss){
      if (kk != jj) {
        # Fir the larger model
        frm <- as.formula(paste0("out ~ (", paste(names(dts.temp)[- 1], 
                                                  collapse= " + "), 
                                 ") * Rgrp"))
        dts.temp$Rgrp <- prep$miss.ind[, prep$locmis[kk]]
        fit.large <- lm(frm, dts.temp)
        # Compare to the larger model
        pval[kk, jj]  <- tryCatch({anova(fit.small, fit.large)[[6]][2]},
                                  error = function (e) NA)
      }
    }
  }
  # Apply the BF correction for multiple testing
  alpha <- alpha / (max(prep$nvar.miss - 1, 1) * prep$nvar.miss)
  reject <- any(pval < alpha, na.rm = TRUE)
  which.reject <- apply(pval, 1, min, na.rm = TRUE) < alpha
  diagMAAR <- list(reject = reject, results = pval, 
                   which.reject = which.reject, method = "ccm")
  class(diagMAAR) <- "diagMAAR"
  return(diagMAAR)
}