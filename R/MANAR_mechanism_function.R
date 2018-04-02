#' MAAR missing data mechanism
#' 
#' \code{MAAR_mechanism} takes a matrix and deletes some of the entries
#' using a missing always not at random mechanism.
#' 
#' \code{MAAR_mechanism} is part of the data generating functions used when
#' conducting simulation studies involving missing data. The probability that 
#' a variable is missing is depends on all variables. 
#' For example, assume we have 3 variables, \eqn{Y_1} and \eqn{Y_2} can have
#' missing values and \eqn{Y_3} is always fully observed. Let \eqn{R_1} and 
#' \eqn{R_2} be the response indicators for \eqn{Y_1} and \eqn{Y_2} respectively
#'  Then the probability that the \eqn{i}-th variable is missing is
#' \deqn{p(R_{i,1}=1|Y_{i,\cdot}) = \text{logit}^{-1}(\alpha_1+\beta_1Y_{i,3}).}
#' \deqn{p(R_{i,2}=1|Y_{i,\cdot},R_{i,1}) = \text{logit}^{-1}(\alpha_2+
#' \gamma_1Y_{i,1}+\beta_2Y_{i,3}).}
#' The vector \eqn{\beta} is specified using the \code{prob.coef} parameter,
#' and \eqn{\alpha} is selected to ensure that the proportion of missing values
#' in each variable is equal to the \code{miss.coef} parameter. 
#' 
#' @param samples A numeric matrix with no missing values.
#' @param miss.coef: A numeric value or a vector of length equal to miss.nvar.
#'              The missingness coefficient(s) determines the proportion of
#'              missing values for each variable.
#' @param miss.nvar An integer specifying the number of variables that will have
#' missing values. miss.nvar <= ncol(samples).
#' @param miss.var An integer vector of length equal to miss.nvar specifying 
#' which variables will have missing values. If left blank a random
#' sample will be taken from the columns of samples will be taken.
#' @param prob.coef A numeric matrix that represents the regression coefficients
#' that will be used to generate the missing data pattern. The
#' nrow(miss.coef) == miss.var and ncol(miss.coef) == ncol(samples). 
#' @param self.dep A logical values indicating if the probability of 
#' variable j missing depends on the values that variable j took. 
#' @param seed A numeric value used to set a seed, leave blank for NULL.
#' @return samples.obs A data matrix of dimensions equal to samples with some
#' missing values.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.mvn <- MANAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 1, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7,- 0.2), 1, 4))
#' @family data generating functions
MANAR_mechanism <- function(
    samples, miss.coef, miss.nvar, miss.var = NULL, prob.coef, 
    self.dep = FALSE, seed = NULL) {
  if (!is.null(seed)) {
    set.seed(seed)
  }
  N <- nrow(samples)
  if (is.null(miss.var)) {
    # If miss.var is not specified randomly sample which columns will have
    # missing values
    miss.var <- sample(ncol(samples), miss.nvar)
  }
  # If only one missingness coefficient is supplied then turn it into a vector
  # for ease of use.
  if (length(miss.coef) != miss.nvar) {
    miss.coef <- rep(miss.coef[1], miss.nvar)
  }
  samples.obs <- samples
  for (kk in 1:miss.nvar) {
    # Compute the probability of missingness
    if (self.dep) {
      # If self.dep is true then the probability of variable j missing depends
      # the values variables j took.
      Xbeta <- samples %*% prob.coef[kk, ]
    } else {
      # Otherwise the probability of the jth variable missing only depends on 
      # the -j variables. 
      Xbeta <- samples[, - miss.var[kk]] %*% prob.coef[kk, ]
    }
    # Use uniroot to find the values of the intercept that makes the mean
    # number of missing values equal to miss.coef
    beta0 <- uniroot(
        logistic_prob, c( - 1000, 1000), Xbeta = Xbeta, mu = miss.coef[kk],
        extendInt = "yes")$root
    # Compute the probability of missing
    p.miss <- logistic_prob(beta0, Xbeta)
    R <- sapply(1:N, function(jj) sample(c(NA, 1), size = 1,
                                         prob = c(p.miss[jj], 1 - p.miss[jj])))
    
    samples.obs[, miss.var[kk]] <- samples[, miss.var[kk]] * R
  }
  return(samples.obs)
}
