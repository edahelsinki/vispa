#' Create null distribution and test statistic for testing patterns on a 
#' banana plot
#' 
#' Create a function f(X) such that: 
#' Given a (n,m) matrix `X`, f(X) evaluates `fun` on `n_lines` subsets of `X`. 
#' The subsets of `X` are obtained by drawing `n_lines` lines on `X`. 
#' The lines are drawn by sampling a starting and an ending point randomly on 
#' the given `grid`, which contains row and column indices of `X` allowed in
#' sampling. 
#' Additional constraints are imposed on grid depending on `line_type`
#' @param smear_df data frame whose first column is a timestamp and the rest 
#' are particle concentrations. Contains multiple days of measurements.
#' @param i_null row indices of smear_df that are null data 
#' @param test_stat_params 
#' list(fun, line_type, n_lines), where
#' * fun is a function, see [make_line_test_stat()]
#' * line_type is one of "rising", "horizontal", "random"
#' * n_lines is the number of lines to sample on the grid
#' @return list(test_stat, sample_null)
#' @export
#'
#' @examples
make_banana_test_params <- function(smear_df, i_null, test_stat_params) {
  num_10min_measurements_in_24hours = 24*60/10 # Note: 23% of days are not complete
  sample_days_smear = make_sample_days(X = smear_df[i_null, ], 
                                       day_length = num_10min_measurements_in_24hours, 
                                       return_matrix = TRUE)
  sample_lines = make_sample_line_segments(smear_df[1:num_10min_measurements_in_24hours,], 
                                           test_stat_params$line_type)
  lines_to_test = sample_lines(test_stat_params$n_lines)
  test_stat = make_test_stat_on_lines(test_stat_params$fun, lines_to_test)
  
  null_samples = split(smear_df[i_null, ], strftime(smear_df[i_null, ]$timestamp, "%Y-%m-%d"))
  # Only 70 days are full days of 144 10-minute measurements, so fill in the missing measurements with NAs.
  i_complete_days = sapply(null_samples, nrow) == 144
  times_in_complete_day = null_samples[[which(i_complete_days)[1]]]$time
  null_samples = lapply(null_samples, function(df) fill_missing_times(df, times_in_complete_day))
  null_samples = null_samples[sapply(null_samples, nrow) == num_10min_measurements_in_24hours & 
                                sapply(null_samples, function(x) !all(is.na(x)))]
  null_samples = lapply(null_samples, function(df) as.matrix(df[,-1]))
  list(sample_null = sample_days_smear, 
       null_samples = null_samples, 
       test_stat = test_stat, 
       test_stat_params = c(test_stat_params, list(lines_to_test = lines_to_test, 
                                                   sample_lines = sample_lines)))
}

#' Make a test statistic for testing banana patterns on a banana plot
#' 
#' Returns a function f(X) that evaluates a function `fun` on X[i,j], where 
#' the indices i and j are determined by drawing multiple lines 
#' `lines_to_test` on the (n, m) grid defined by the rows and columns of X.
#' 
#' @param fun function fun(x) that returns a single numeric value, 
#' given a numeric vector x
#' @param lines_to_test (n_lines, 4) matrix whose columns are the starting and 
#' end indices (x0,y0,x1,y1) for lines drawn on a (n,m) grid. They must follow the 
#' constraints: 1<=x0,x1<=m and 1<=y0,y1<=n.
#'
#' @return f(X) where X is a (n,m) matrix
#' @export
#'
#' @examples
make_test_stat_on_lines = function(fun, lines_to_test) {
  function(X) apply(lines_to_test, 1, 
                    function(line) fun(get_matrix_values_on_line(X, line)))
}

#' Make function for sampling lines from a template
#'
#' @param df data frame, first column is timestamp
#' @param line_type one of "rising", "horizontal", "random"
#'
#' @return f(n_lines) that returns (n_lines, 4) matrix
#' @seealso [sample_line_segments()]
#' @export
#'
#' @examples
make_sample_line_segments = function(df, line_type="rising", grid = NULL) {
  if (is.null(grid))
    grid = make_grid_indices(dim(df[, -1]))
  j_less_than_25nm <- which(particle_colname_to_nanometers(colnames(df[, -1])) < 25)
  one_hour = 6 # one hour grid spacing for 10-minute measurements
  
  if (line_type=="rising") {
    i_start = grid$x[1:(length(grid$x)-one_hour-1)]
    j_start = intersect(grid$y, j_less_than_25nm)
    make_i_end = function(p0) grid$x[grid$x - p0[1] > one_hour] # x1-x0>1h
    make_j_end = function(p0) grid$y[grid$y > p0[2]] # y1>y0
  } else if (line_type=="horiz") {
    i_start = grid$x[-length(grid$x)]
    j_start = grid$y
    make_i_end = function(p0) grid$x[grid$x > p0[1]] # x1>x0
    make_j_end = function(p0) grid$y[grid$y == p0[2]] # y1==y0
  } else if (line_type=="random") {
    i_start = grid$x
    j_start = grid$y
    make_i_end = function(p0) grid$x
    make_j_end = function(p0) grid$y
  } 
  
  function(n_lines) sample_line_segments(i_start, j_start, n_lines, 
                                         make_i_end, make_j_end) 
}

#' Sample line segments on a rectangular grid
#'
#' @param i,j Allowed row and column indices of a 2D grid from which to sample starting 
#' points
#' @param make_i_end,make_j_end functions f(p0) that return the allowed indices 
#' of a 2D grid for sampling the lines' ending points, given the starting point 
#' p0 = c(x0, y0). 
#' @param size number of lines to sample
#'
#' @return (size, 4) matrix of (x0, y0, x1, y1) coordinates for starting points 
#' p0=(x0,y0) and ending points p1=(x1,y1) that create a line segment, where 
#' p0 != p1
#' @export
#'
#' @examples
#' sample_line_segments(1:5, 1:4, 4)
#' 
#' # Visualisation of what is happening
#' set.seed(42)
#' X = matrix(rnorm(100*100), 100, 100)
#' L = sample_line_segments(1:nrow(X), 1:nrow(X), 10)
#' image(1:nrow(X), 1:ncol(X), X)
#' for (i in 1:nrow(L)) 
#'   segments(L[i, 1], L[i, 2], L[i, 3], L[i, 4])
sample_line_segments <- function(i, j, size = 1, make_i_end = NULL, make_j_end = NULL) {
  if (is.null(make_i_end))
    make_i_end = function(p0) i
  if (is.null(make_j_end))
    make_j_end = function(p0) j
  points_start = sample_grid_points(i, j, size)
  points_end = points_start
  for (k in 1:size) {
    while (identical(points_start[k,], points_end[k,]))
      points_end[k,] = sample_grid_points(make_i_end(points_start[k,]), 
                                          make_j_end(points_start[k,]), 
                                          1)
  }
  lines = cbind(points_start, points_end)
  colnames(lines) = c("x0", "y0", "x1", "y1")
  lines
}

#' Sample a point from a rectangular grid
#' @param i,j Row and column indices of a 2D grid.
#' @return (size, 2) matrix of x,y grid coordinates for sampled points
#' @export
#' sample_grid_points(1:5, 1:4, 4)
sample_grid_points <- function(i, j, size=1) {
  stopifnot(length(i) > 0, length(j) > 0)
  # Why check length(i) == 1: 
  # If i or j are a single integer, then the default behavior of sample(i) 
  # means sample(1:i).
  # But i,j denote the grid indices from which it is allowed to sample, so 
  # sample(i, size) should return rep(i, size). 
  # An example of where this is a problem is for example when there is a 
  # constraint of sampling only horizontal lines, where j=y0=y1.
  x = if (length(i) == 1) rep(i, size) else sample(i, size, replace=T)
  y = if (length(j) == 1) rep(j, size) else sample(j, size, replace=T)
  matrix(c(x, y), ncol=2, dimnames=list(NULL, c("x","y")))
}

#' Create a rectangular grid
#'
#' @param sizes integer vector of grid sizes
#' @param spacing integer vector of grid spacing. Spacing for horizontal and vertical direction. 
#' For example, c(3,2) keep every 3rd point in x and every 2nd point in y. 
#' @return list(x, y) where x and y are integer vectors, if length(sizes)=2. 
#' Else list of same length as sizes
#' @examples
#' # Creates indices for 2D grid. Equivalent to list(x=1:10, y=1:15)
#' make_grid_indices(c(10, 15))
#' # Creates indices for 2D grid where x=1:10 spaced by 2 and y=1:15 spaced by 3
#' make_grid_indices(c(10, 15), c(2,3))
#' # Creates 3D grid
#' make_grid_indices(c(10, 15, 3), c(2,3,1))
make_grid_indices <- function(sizes, spacing = c(1, 1)) {
  stopifnot(length(sizes) == length(spacing))
  grid = vector("list", length(sizes))
  if (length(sizes) == 2) {
    names(grid) = c("x", "y")
  } else {
    names(grid) = paste0("v", 1:length(sizes))
  }
  for (i in seq_along(sizes)) 
    grid[[i]] = seq(1, sizes[i], spacing[i])
  grid
}

#' Get matrix values on line segment
#' Draw a straight line on a matrix, round the line on the discrete matrix grid, 
#' and take the values of the matrix on those points.
#' This is a simplified "vector-to-raster" conversion.
#' @param mat (n,m) matrix
#' @param line vector of indices of mat c(x0, y0, x1, y1) 
#' where 1<x0,x1< n and 1<y0,y1 < m
#' @export
get_matrix_values_on_line = function(mat, line) {
  stopifnot(is.matrix(mat))
  x0 = line[1]
  y0 = line[2]
  x1 = line[3]
  y1 = line[4]
  i = x0:x1
  j = if (x0 == x1) y0:y1 else round((i - x0) * (y1 - y0) / (x1 - x0) + y0)
  diag(mat[i, j, drop=F])
}

#' Fill in missing times of the day with NAs
#'
#'
#' @param df data frame with column timestamp
#' @param times_in_complete_day vector of times in a complete day in format "%H:%M:%S"
#'
#' @return df, with filled in `times_in_complete_day`
#' @export
#'
#' @examples
fill_missing_times = function(df, times_in_complete_day) {
  Mode <- function(x) {
    ux <- unique(x)
    ux[which.max(tabulate(match(x, ux)))]
  }
  date = Mode(strftime(df$timestamp, "%Y-%m-%d"))
  times = strftime(df$timestamp, "%H:%M:%S")
  if (length(times) == length(times_in_complete_day))
    return(df)
  missing_times = setdiff(times_in_complete_day, times)
  df_supplement = data.frame(matrix(nrow=length(missing_times), 
                                    ncol=ncol(df), 
                                    dimnames=list(NULL, colnames(df))))
  df_supplement$timestamp = strptime(paste(date, missing_times), format="%Y-%m-%d %H:%M:%S")
  df_out = rbind(df, df_supplement)
  df_out[order(df_out$timestamp), ]
}

#' Convert column names in SMEAR data into numeric, e.g. d282e1 into 2.82
#'
#' @param smear_names colnames(banana_data) minus timestamp
#'
#' @return numeric vector
particle_colname_to_nanometers <- function(smear_names) {
  s = gsub("d", "", smear_names, fixed=T)
  as.numeric(s) / 1e3
}
