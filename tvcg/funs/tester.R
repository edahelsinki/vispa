#' Perform test
#' Compute `test_stat()` on `X_obs` and on surrogates generated from `sample_from_null()`.
#'
#' @param X_obs             (n,m) matrix of observed data
#' @param n_resample        number of samples from null
#' @param test_stat         function that accepts (n,m) matrix and outputs a vector of k test statistics
#' @param sample_from_null  either: 
#' a function that outputs a data sample from the null as a (n,m) matrix, or
#' a function that receives `n_resample` as a first argument and outputs `n_resample` test statistic samples from the null as a (k, n_resample) matrix
#' @param samples_test_statistics      if TRUE, then sample_from_null has `n_resample` as first argument and returns test statistic samples instead of data samples
#' @param ...               additional parameters for test_stat
#'
#' @return
#' @export
tester <- function(X_obs, 
                   n_resample = 1000, 
                   test_stat, 
                   sample_from_null, 
                   samples_test_statistics=F, 
                   ...) {
  test_stat0 <- function(x) test_stat(x, ...)
  t0 <- test_stat0(X_obs)
  if (samples_test_statistics) {
    t_null = sample_from_null(n_resample)
  } else {
    t_null <- replicate(n_resample, test_stat0(sample_from_null()))
    if (is.list(t_null)) 
      stop("Replicate output is a list. Check output of test_stat (e.g. NULL causes replicate to return a list).")
  }
  list(t0 = t0,
       t_null = t_null,
       params = list(test_stat=test_stat, 
                     sample_from_null=sample_from_null, 
                     test_stat_params=list(...)))
}

#' Tester for given null samples
#'
#' @param X_obs data frame or matrix
#' @param null_samples list null samples
#' @param test_stat function to be applied to X_obs and to elements of null_samples
#'
#' @return
#' @export
#'
#' @examples
tester2 <- function(X_obs, null_samples, test_stat) {
  t0 = test_stat(X_obs)
  t_null = vapply(null_samples, test_stat, t0)
  list(t0 = t0,
       t_null = t_null, 
       params = list(test_stat=test_stat, 
                     null_samples=null_samples))
}

