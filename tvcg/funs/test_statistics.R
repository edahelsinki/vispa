library(sp)
library(scagnostics)

#' Test statistic: scagnostics
#'
#' @param data    nX2 matrix
#'
#' @return        1X9 vector of scagnostics for `data` in view defined by `pvec`
#' @export
scag <- function(xy) {
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

diff_in_interval <- function(x, t1, t2) x[t2] - x[t1]

diff_in_all_2h_intervals <- function(x) apply_to_intervals(as.matrix(x),
                                                           f = diff_in_interval,
                                                           L = 2,
                                                           returnNamed = T)

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

#' Periodicity (spectrum)
#'
#' @param x vector
#'
#' @return numeric vector of frequency magnitudes
#' @export
#'
#' @examples
periodicity = function(x) {
  x_F = abs(fft(x)) 
  x_F[1:floor(length(x_F)/2)]
}

#' Names of frequencies for `periodicity()`
#'
#' @param x numeric vector
#' @param f_s sampling rate
#'
#' @return numeric vector of frequencies
#' @export
#'
#' @examples
periodicity_names = function(x, f_s=1) {
  f=0:(length(x)-1)/length(x)*f_s 
  f[1:floor(length(f)/2)] 
} 

#' Histogram bins
#'
#' @param x 
#' @param breaks 
#' @param norm 
#'
#' @return
#' @export
#'
#' @examples
hist_test_stat = function(x, breaks, norm=TRUE) {
  y=tabulate(cut(x, breaks), length(breaks)-1)
  if (norm) y/sum(y) else y
}

#' Names for histogram bins
#'
#' @param x 
#' @param breaks 
#'
#' @return
#' @export
#'
#' @examples
hist_test_stat_names = function(x, breaks) {
  levels(cut(x, breaks))
}
