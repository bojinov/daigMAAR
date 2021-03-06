% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/MAAR_mechanism_function.R
\name{MAAR_mechanism}
\alias{MAAR_mechanism}
\title{MAAR missing data mechanism}
\usage{
MAAR_mechanism(samples, miss.coef, miss.nvar, miss.var = NULL, prob.coef,
  seed = NULL)
}
\arguments{
\item{samples}{A numeric matrix with no missing values.}

\item{miss.nvar}{An integer specifying the number of variables that will have
missing values. miss.nvar < ncol(samples).}

\item{miss.var}{An integer vector of length equal to miss.nvar specifying 
which variables will have missing values. If left blank a random
sample will be taken from the columns of samples will be taken.}

\item{prob.coef}{A numeric matrix that represents the regression coefficients
that will be used to generate the missing data pattern. The
nrow(miss.coef) == miss.var and ncol(miss.coef) == ncol(samples) - miss.var.}

\item{seed}{A numeric value used to set a seed, leave blank for NULL.}

\item{miss.coef:}{A numeric value or a vector of length equal to miss.nvar.
The missingness coefficient(s) determines the proportion of missing values
for each variable.}
}
\value{
samples.obs A data matrix of dimensions equal to samples with some
missing values.
}
\description{
\code{MAAR_mechanism} takes a matrix and deletes some of the entries
using a missing always at random mechanism.
}
\details{
\code{MAAR_mechanism} is part of the data generating functions used when
conducting simulation studies involving missing data. The probability that 
a variable is missing is only a function of the fully observed variables. 
For example, assume we have 3 variables, \eqn{Y_1} and \eqn{Y_2} can have
missing values and \eqn{Y_3} is always fully observed. Let \eqn{R_1} and 
\eqn{R_2} be the response indicators for \eqn{Y_1} and \eqn{Y_2} respectively
 Then the probability that the \eqn{i}-th variable is missing is
\deqn{p(R_{i,1}=1|Y_{i,\cdot}) = \text{logit}^{-1}(\alpha_1+\beta_1Y_{i,3}).}
\deqn{p(R_{i,2}=1|Y_{i,\cdot}) = \text{logit}^{-1}(\alpha_2+\beta_2Y_{i,3}).}
The vector \eqn{\beta} is specified using the \code{prob.coef} parameter,
and \eqn{\alpha} is selected to ensure that the proportion of missing values
in each variable is equal to the \code{miss.coef} parameter.
}
\examples{
# Generate 100 iid samples from a MVN with correlation equal to 0.3
samples.mvn <- sample_mvn(5, 0.3, 100)
obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
                           miss.nvar = 1, miss.var = 4, 
                           prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
}
\seealso{
Other data generating functions: \code{\link{MAAR_mechanism_noindep}},
  \code{\link{MANAR_mechanism}},
  \code{\link{logistic_prob}}, \code{\link{sample_mvn}}
}
