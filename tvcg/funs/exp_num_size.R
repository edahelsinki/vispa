source("conttab.R")
source("null_distributions.R")
source("pvalue.R")
#source("conttab_power_experiment.R")
source("conttab_toydata.R")

#' Do contingency table experiment for various sizes
#'
#' @param n_grid 
#' @param m_grid 
#' @param num_tile_grid 
#' @param size_tile_grid 
#' @param verbose 
#' @param n_resample 
#'
#' @return data frame of every combination of n_grid, m_grid, num_ile_grid, size_tile_grid, 
#' and an experiment column where experiment = list(X_obs, tiles_planted, pv)
#' @export
#'
#' @examples
#' set.seed(2020)
#' exp1 = do_experiment_num_size(1000, 6, c(10, 100), c(10, 50), 1000, TRUE)
do_experiment_num_size = function(n_grid, 
                                  m_grid, 
                                  num_tile_grid,
                                  size_tile_grid, 
                                  n_resample = 1000, 
                                  verbose=FALSE) {
  exp_df = expand.grid(n = n_grid, 
                       m = m_grid, 
                       size_tile = size_tile_grid, 
                       num_tile = num_tile_grid, 
                       experiment = list(NULL), 
                       KEEP.OUT.ATTRS = FALSE)
  exp_df = exp_df[exp_df$n > exp_df$size_tile & 
                    exp_df$n*exp_df$m > 2*exp_df$num_tile*exp_df$size_tile, ]
  for (i in 1:nrow(exp_df)) {
    n = exp_df$n[i]
    m = exp_df$m[i]
    num_tile = exp_df$num_tile[i]
    size_tile = exp_df$size_tile[i]
    if (verbose)
      cat(paste0(i, "/", nrow(exp_df), 
                 "\t n=", n, 
                 "\t m=", m, 
                 "\t num_tile=", num_tile, 
                 "\t size_tile=", size_tile, "\n"))
    tiles_planted = make_non_overlapping_tiles(n, m, size_tile, num_tile)
    X_obs = make_data_with_patterns_two_column(n, m, tiles_planted)
    pv = test_many_conttab(X_obs, tiles_planted, n_resample)
    exp_df$experiment[[i]] = list(X_obs=X_obs, tiles_planted=tiles_planted, pv=pv)
  }
  exp_df
}

#' Test multiple contingency tables
#' 
#' Compute p-values for the cells of `table(X_obs[R, C])`, 
#' for tiles (R,C) in tiles_to_test. 
#'
#' @param X_obs data frame or matrix
#' @param tiles_to_test matrix with columns r0, r1, C1, C2
#' @param n_resample number of null samples for p-value computation
#'
#' @return pvalues for every pattern in tiles_to_test
#' @export
#'
#' @examples
#' 
test_many_conttab <- function(X_obs, tiles_to_test, n_resample) {
  tiling_H0 = tiling(n=nrow(X_obs), m=ncol(X_obs))
  n_pattern = nrow(tiles_to_test)
  pv = matrix(nrow=n_pattern, ncol=4)
  for (i in 1:n_pattern) {
    tile = tiles_to_test[i, ]
    R = tile['r0']:tile['r1']
    C = c(tile['C1'], tile['C2'])
    test = tester_conttab(X_obs, tiling_H0, n_resample, R, C)
    pv[i,] = test_to_pvalue(test, use_minp=TRUE)
  }
  pv
}

#' Enumerate non-overlapping tiles on a matrix
#' 
#' Enumerate num_tile non-overlapping tiles (R,C) of size size_tile 
#' where R=r0:r1, C=c(C1, C2). 
#' Starts from R=1:size_tile, C=c(1,2).
#'
#' @param n,m 
#' @param size_tile,num_tile 
#'
#' @return integer matrix with "r0", "r1", "C1", "C2"
#' @export
#'
#' @examples
#' make_non_overlapping_tiles(10, 4, 2, 8)
make_non_overlapping_tiles <- function(n, m, size_tile, num_tile) {
  tiles = matrix(0L, nrow=num_tile, ncol=4, dimnames=list(NULL, c("r0", "r1", "C1", "C2")))
  size_tile = as.integer(size_tile)
  r0 = 1L
  r1 = size_tile
  C1 = 1L
  C2 = 2L
  for (i in 1:nrow(tiles)) {
    tiles[i, ] = c(r0, r1, C1, C2)
    C1 = C1 + 2L
    C2 = C2 + 2L
    if (C2 > m) {
      C1 = 1L
      C2 = 2L
      r0 = r1 + 1L
      r1 = r0 + size_tile - 1L
    }
    if (r1 > n) {
      warning(paste0("Can't sample ", num_tile, " non-overlapping tiles of size ", size_tile, 
                     " from a (", n, ",", m, ") matrix.", "\n", 
                     "Using ", i, " tiles instead."))
      tiles = tiles[1:i, ]
      break 
    }
  }
  tiles
}
