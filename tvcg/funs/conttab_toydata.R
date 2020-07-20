source("null_distributions.R")

#' Make toy data
#' 
#' Make `n*m` data frame with binary factor columns. The columns are fully 
#' correlated in rows and column indices given in `patterns`.
#'
#' @param n,m
#' @param patterns list of patterns, where each pattern is a list(R, C), where 
#' R are row indices and C are column indices
#'
#' @return data frame with binary factor columns
#' @export
#'
#' @examples
#' # Make data where the first four rows are correlated: 
#' make_data_with_patterns(10, 2, list(list(R=1:4, C=1:2)))
#' make_data_with_patterns(20, 3, list(list(R=1:4, C=1:2), list(R=9:15, C=2:3)))
make_data_with_patterns = function(n, m, patterns) {
  X = data.frame(matrix(nrow=n, ncol=m))
  for (j in 1:m)
    X[,j] = factor(rep(c(0,1), length.out=n), levels=c(0,1), ordered = FALSE)
  tiling_H1 = tiling(n, m)
  if (max(unlist(lapply(patterns, function(x) x$R))) > n || 
      max(unlist(lapply(patterns, function(x) x$C))) > m)
    stop("Some R or C in patterns is larger than n or m.")
  for (i in seq_along(patterns))
    tiling_H1$addtile(R = patterns[[i]]$R, C = patterns[[i]]$C)
  X_obs = tiling_H1$permutedata(X)
  attr(X_obs, "patterns") = patterns
  X_obs
}

#' Make toy data
#' 
#' Wrapper function for [make_data_with_patterns()].
#'
#' @inheritParams n,m 
#' @param pattern_df data frame with columns r0,r1,C1,C2 
#' @examples 
#' make_data_with_patterns_two_column(10, 2, data.frame(r0=1, r1=4, C1=1, C2=2))
make_data_with_patterns_two_column = function(n, m, patterns_df) {
  patterns = apply(patterns_df, 1, 
                   function(p) list(R=p["r0"]:p["r1"], C=unname(c(p["C1"], p["C2"]))))
  make_data_with_patterns(n, m, patterns)
}

#' Enumerate all possible 2-column tiles
#'
#' @param n,m number of rows and columns in tiling
#' @param size_patterns if non-zero, returns only patterns of height `size_patterns`
#'
#' @return data frame with r0, r1, C1, C2, where r0,r1 are standing and ending 
#' row indices, and C1,C2 are column indices, such that: 
#' `tile = (R, C) = (r0:r1, (C1,C2))`
#' @export
#'
#' @examples
enumerate_two_column_tiles <- function(n, m, size_patterns=0) {
  pattern_space = expand.grid(r0=1:(n-1), 
                              r1=2:n, 
                              C1=1:m, 
                              C2=1:m, 
                              KEEP.OUT.ATTRS = FALSE)
  pattern_space = pattern_space[pattern_space$C1 < pattern_space$C2 & 
                                  pattern_space$r0 < pattern_space$r1, ]
  if (size_patterns != 0) 
    pattern_space = pattern_space[(pattern_space$r1 - pattern_space$r0 + 1) == size_patterns, ]
  pattern_space = pattern_space[order(pattern_space$r0, 
                                      pattern_space$C1,
                                      pattern_space$C2), ]
  
  rownames(pattern_space) = NULL
  pattern_space
}
