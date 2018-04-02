#' Prepare for diagnostic
#' 
#' \code{prep} takes a data frame, or matrix and prepares it to be analyzed
#'  by diagMAAR. 
#' 
#' \code{prep} is part of the data preparation functions used to pre-process
#' data sets so that missing data diagnostic tests can easily be applied.
#' Currently there are three diagnostic test implemented: A comparison of 
#' conditional means (ccm), directly testing a postulated missingness 
#' mechanism (dtmm), and Gaussian copula (cop). 
#' 
#' @param dt A numeric matrix, data frame or data table.
#' @param daigtest A string indicating which diagnostic daigtest to use. The
#' options are: cop, dtmm, ppp. See help for prep.daigtest.
#' @return prep A prep S3  object that contains the original data,
#' the desired daigtest to use, the sample size, the number of missing
#' variables, the location of the missing variables, the number of variables and
#' a data set that has the missingness indicators appended.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 4, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
#' dt.prep <- prep(obs.mvn, "cop")
#' print(dt.prep)
#' dt.prep <- prep(obs.mvn, "dtmm")
#' print(dt.prep)
#' dt.prep <- prep(obs.mvn, "ccm")
#' print(dt.prep)
#' @family prep
prep <- function(dt, daigtest) {
  if (daigtest == "cop") {
    prep <- prep.cop(dt)
  } else if (daigtest == "dtmm") {
    prep <- prep.dtmm(dt)
  } else if (daigtest == "ccm") {
    prep <- prep.ccm(dt) 
  } else {
    cat("Invalid daigtest specified. Please use cop, dtmm, or ccm only")
    prep <- NULL
  }
  return(prep)
}

#' Prepare for copula diagnostics
#' 
#' \code{prep.cop} takes a data frame, or matrix and prepares it to be analyzed
#'  by diagMAAR.cop.
#' 
#' \code{prep.cop} is part of the data preparation functions used to pre-process
#' data sets so that diagnostic tests for the MAAR assumption can be ran. 
#' 
#' @param dt A numeric matrix, data frame or data table.
#' @return prep A prep S3  object that contains the original data,
#' the desired daigtest to use, the sample size, the number of missing
#' variables, the location of the missing variables, the number of variables and
#' a data set that has the missingness indicators appended.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 4, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
#' dt.prep <- prep.cop(obs.mvn)
#' print(dt.prep)
#' @family prep
prep.cop <- function(dt) {
  # Extract the missing data indicators.
  dt <- as.matrix(dt)
  N <- dim(dt)[1]
  nvar <- dim(dt)[2]
  # Generate the prep
  prep <- list(data.original = dt, daigtest = "cop", N = N, nvar = nvar)
  miss.ind <- is.na(dt) + 0
  prep$miss.ind <- miss.ind
  prep$locmis <- which(colSums(miss.ind) > 0)
  prep$nvar.miss <- length(prep$locmis)
  # Get which variable have missing values.
  full.obs.cols <- which(colSums(miss.ind) == 0)
  # Combine the data sat and the missingness indicators.
  dt <- cbind(dt[, - full.obs.cols], miss.ind[, - full.obs.cols], 
              dt[, full.obs.cols])
  prep$dt <- dt
  class(prep) <- "prep"
  return(prep)
}

#' Prepare for directly testing a postulated missingness mechanism
#' 
#' \code{prep.dtmm} takes a data frame, or matrix and prepares it to be analyzed
#'  by diagMAAR.dtmm.
#' 
#' \code{prep.dtmm} is part of the data preparation functions used to pre-process
#' data sets so that diagnostic tests for the MAAR assumption can be ran. 
#' 
#' @param dt A numeric matrix, data frame or data table.
#' @return prep A prep S3  object that contains the original data,
#' the desired daigtest to use, the sample size, the number of missing
#' variables, the location of the missing variables, the number of variables and
#' a data set that has the missingness indicators appended.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 4, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
#' dt.prep <- prep.dtmm(obs.mvn)
#' print(dt.prep)
#' @family prep
prep.dtmm <- function(dt) {
  # Extract the missing data indicators.
  dt <- as.matrix(dt)
  N <- dim(dt)[1]
  nvar <- dim(dt)[2]
  # Generate the prep
  prep <- list(data.original = dt, daigtest = "dtmm", N = N, nvar = nvar)
  miss.ind <- is.na(dt) + 0
  prep$miss.ind <- miss.ind
  prep$locmis <- which(colSums(miss.ind) > 0)
  prep$nvar.miss <- length(prep$locmis)
  # Get which variable have missing values.
  full.obs.cols <- which(colSums(miss.ind) == 0)
  # Combine the data sat and the missingness indicators.
  dt <- cbind(dt, miss.ind[, - full.obs.cols])
  colnames(dt) <- c(sapply(1:nvar, 
                           function(ll) paste0("Y", as.character(ll))),
                    sapply(1:prep$nvar.miss , 
                           function(ll) paste0("R", as.character(ll))))
  prep$dt <- dt
  class(prep) <- "prep"
  return(prep)
}
#' Prepare for comparison of conditional means
#' 
#' \code{prep.ccm} takes a data frame, or matrix and prepares it to be analyzed
#'  by diagMAAR.ccm.
#' 
#' \code{prep.ccm} is part of the data preparation functions used to pre-process
#' data sets so that diagnostic tests for the MAAR assumption can be ran. 
#' 
#' @param dt A numeric matrix, data frame or data table.
#' @return prep A prep S3  object that contains the original data,
#' the desired daigtest to use, the sample size, the number of missing
#' variables, the location of the missing variables, the number of variables and
#' a data set that has the missingness indicators appended.
#' @export
#' @examples 
#' # Generate 100 iid samples from a MVN with correlation equal to 0.3
#' samples.mvn <- sample_mvn(5, 0.3, 100)
#' # Take the Gaussian data and and delete some values from the fourth row.
#' obs.nvmn <- MAAR_mechanism(samples = samples.mvn, miss.coef = 0.2, 
#'                            miss.nvar = 1, miss.var = 4, 
#'                            prob.coef = matrix(c(-1, 0.5, 0.7, - 0.2), 1, 4))
#' dt.prep <- prep.ccm(obs.mvn)
#' print(dt.prep)
#' @family prep
prep.ccm <- function(dt) {
  # Extract the missing data indicators.
  dt <- as.matrix(dt)
  N <- dim(dt)[1]
  nvar <- dim(dt)[2]
  # Generate the prep
  prep <- list(data.original = dt, daigtest = "ccm", N = N, nvar = nvar)
  miss.ind <- is.na(dt) + 0
  prep$miss.ind <- miss.ind
  prep$locmis <- which(colSums(miss.ind) > 0)
  prep$nvar.miss <- length(prep$locmis)
  # Get which variable have missing values.
  full.obs.cols <- which(colSums(miss.ind) == 0)
  # Combine the data sat and the missingness indicators.
  prep$dt <- dt
  class(prep) <- "prep"
  return(prep)
}
#' @export
print.prep <- function(prep, ...) {
  cat("daigtest:")
  print(prep$daigtest)
  cat("\nRows:")
  print(prep$N)
  cat("\nNumber of Variables:")
  print(prep$nvar)
  cat("\nNumber of Variables With Missing Values:")
  print(prep$nvar.miss)
  cat("\nTable of Missing Variables:")
  print(colSums(prep$miss.ind))
}
