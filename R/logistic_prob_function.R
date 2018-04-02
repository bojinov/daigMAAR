#' Logistic Probability
#' 
#' This function computes \eqn{logit^{-1}(\beta_0 + X\beta)} if is.null(mu), 
#' otherwise it computes \eqn{mean(logit^{-1}(\beta_0 + X\beta) - \mu)}.
#' 
#' This function is used when generating missing data from a complete data set.
#' 
#' 
#' @param beta0: A numeric value.
#' @param Xbeta: A numeric matrix.
#' @param mu: A numeric value, by default it is set to NULL
#' @return p: logit^{-1}(beta0 + Xbeta) if is.null(mu) or 
#' mean(logit^{-1}(beta0 + Xbeta)) - mu.
#' @family data generating functions
#' @export

logistic_prob <- function(beta0, Xbeta, mu = NULL) {
  if (is.null(mu)) {
    p <- 1 / (1 + exp( - (beta0 + Xbeta)))
  } else {
    p <- mean(1 / (1 + exp( - (beta0 + Xbeta)))) - mu
  }
  return(p)
}
