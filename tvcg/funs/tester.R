library(sp)
library(scagnostics)

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

#### Test statistics ####
 
#' Test statistic: scagnostics
#'
#' @param data    nX2 matrix
#'
#' @return        1X9 vector of scagnostics for `data` in view defined by `pvec`
#' @export
scag <- function(xy) {
  # binostics::scagnostics seems very slow. On average, 380ms for ONE evaluation 
  # of scagnostics on the german data (400 rows).
  # scagnostics::scagnostics on the other hand is 15ms. But it requires rJava.. 
  #as.matrix(binostics::scagnostics(xy))[,]
  #as.matrix(scagnostics::scagnostics(xy))[,]
  scagnostics::scagnostics(xy)
}

#' Test statistic: number of points inside selected region.
#'
#' @param data    nXm matrix
#' @param region  kX2 matrix with x,y coordinates of a polygon selection
#' @param pvec      mX2 matrix with projection vectors (vectors as columns)
#'
#' @return        number of points inside selected region.
#' @export
num_points_in_selection <- function(data, pvec, region) {
  length(points_in_selection(data %*% pvec, region))
}

#' Returns rownames of data points inside a polygon selection. Utilizes sp library. 
#'
#' @param data      nX2 matrix with xy coordinates
#' @param selection kX2 matrix of xy coordinates for polygon vertices 
#'
#' @return
points_in_selection <- function(data, selection) {
  data <- data.frame(data)
  selection <- data.frame(selection)
  idx <- sp::point.in.polygon(point.x = data[ , 1], point.y = data[ , 2], 
                              pol.x = selection[ , 1], pol.y = selection[ , 2])
  if (length(idx) < 2) warning(length(idx), " points in selection.")
  as.integer(rownames(data)[as.logical(idx)])
}

#' Test statistic: Apply function f to intervals [a, b] of x.
#'
#' @param x matrix, array or data.frame nXm
#' @param f function(x, a, b) of interval [a, b]
#' @param L vector of interval sizes. f is applied only on intervals with lengths in L.
#' @param returnNamed if TRUE, return named vector. Name="a_b" for value at [a,b].
#'
#' @return array kXp where k = output size of f, p = # of intervals in x.
#' @examples apply_to_intervals(x = array(1:100), f = function(x, a, b) x[b] - x[a], L = c(4,5)) returns vector of x[b]-x[a] for all intervals [a,b] of length 4 or 5.
#' @export
apply_to_intervals <- function(x, f, L=NULL, returnNamed = F) {
  n <- dim(x)[1]
  x_intervals <- combn(1:n, 2)
  if (!is.null(L)) {
    if (L >= n) stop("L>=n: interval size L cannot be longer than dim(x)[1]")
    colDiffs <- function(x) colSums(x * c(-1,1))
    x_intervals <- x_intervals[ , colDiffs(x_intervals) %in% L, drop=F]
  }
  y <- apply(x_intervals, 2, FUN=function(interval) f(x, interval[1], interval[2]))
  if (returnNamed) names(y) <- apply(x_intervals, 2, function(interval) paste0(interval[1], "-", interval[2]))
  y
}