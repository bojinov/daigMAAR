% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prep_function.R
\name{prep}
\alias{prep}
\title{Prepare for diagnostic}
\usage{
prep(dt, daigtest)
}
\arguments{
\item{dt}{A numeric matrix, data frame or data table.}

\item{daigtest}{A string indicating which diagnostic daigtest to use. The
options are: cop, dtmm, ppp. See help for prep.daigtest.}
}
\value{
prep A prep S3  object that contains the original data,
the desired diagnostic test to use, the sample size, the number of missing
variables, the location of the missing variables, the number of variables and
a data set that has the missingness indicators appended.
}
\description{
\code{prep} takes a data frame, or matrix and prepares it to be analyzed
 by diagMAAR.
}
\details{
\code{prep} is part of the data preparation functions used to pre-process
data sets so that missing data diagnostic tests can easily be applied.
Currently there are three diagnostic test implemented: A comparison of 
conditional means (ccm), directly testing a postulated missingness 
mechanism (dtmm), and Gaussian copula (cop).
}
\examples{
# Generate 100 iid samples from a MVN with correlation equal to 0.3
samples.mvn <- sample_mvn(5, 0.3, 100)
# Take the Gaussian data and and delete some values from the fourth row.
obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
                           miss.nvar = 1, miss.var = 4, 
                           prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
dt.prep <- prep(obs.mvn, "cop")
print(dt.prep)
dt.prep <- prep(obs.mvn, "dtmm")
print(dt.prep)
dt.prep <- prep(obs.mvn, "ccm")
print(dt.prep)
}
\seealso{
Other prep: \code{\link{prep.ccm}}, \code{\link{prep.cop}},
  \code{\link{prep.dtmm}}
}
