source("conttab.R", chdir = TRUE)
source("conttab_toydata.R", chdir = TRUE)
source("null_distributions.R", chdir = TRUE)
source("pvalue.R", chdir = TRUE)

#' Measure signal in overlapping patterns
#' 
#' Since a two-column matrix is used, a tile is parametrised only by its row indices R.
#' 
#' 1. Create observed data `X_obs`: Create a binary (n,2) matrix. Add a tile `R_real=1:size_pattern`, 
#' C=1:2. Permute tiling once to obtain `X_obs`.
#' 2. Perform a test for all possible tiles `R=r0:r1` of maximum size 
#' `max_size_pattern_to_test=r1-r0+1`: compute a contingency table on `X_obs[R,1:2]`, 
#' compute p-values using minP on null data obtained by permuting `X_obs`.
#' 3. Save the minimum p-value in each table, and compute `pct_signal = |R^R_real| / |R|`.
#'
#' @param n number of total rows
#' @param size_pattern size of planted pattern R_real=1:size_pattern
#' @param max_size_pattern_to_test maximum size of tested patterns R=r0:r1, 
#' for all r0<size_pattern,r1-r0<max_size_pattern_to_test
#'
#' @return
#' @export
#'
#' @examples
do_experiment_overlap <- function(n, size_pattern, max_size_pattern_to_test, verbose=FALSE) {
  m = 2
  X_obs = make_data_with_patterns(n, m, list(list(R=1:size_pattern, C=1:2)))
  tiling_H0 = tiling(n, m)
  patterns_to_test = enumerate_two_column_tiles(n, m)
  patterns_to_test = patterns_to_test[patterns_to_test$r0 <= size_pattern & 
                                        patterns_to_test$r1 - patterns_to_test$r0 + 1 <= 
                                        max_size_pattern_to_test, ]
  rownames(patterns_to_test) = NULL
  test_list = vector("list", nrow(patterns_to_test))
  for (i in 1:nrow(patterns_to_test)) {
    if (verbose && (i%%50==0 | i==1))
      cat(paste0("Testing pattern ", i, "/", nrow(patterns_to_test)), "\n")
    R = patterns_to_test$r0[i]:patterns_to_test$r1[i]
    C = c(patterns_to_test$C1[i], patterns_to_test$C2[i])
    test_list[[i]] = tester_conttab(X = X_obs, 
                                    null_tiling = tiling_H0, 
                                    nresample = 1000, 
                                    R = R, 
                                    C = C)
  }
  pv = vapply(test_list, function(test) minP_array(test$t0, test$t_null), 
              matrix(numeric(1), nrow=2,ncol=2))
  pv_min = apply(pv, 3, min)
  pct_signal = apply(patterns_to_test, 1, function(P) {R1=P['r0']:P['r1']; R2=1:size_pattern; length(intersect(R1, R2)) / length(R1)})
  list(result=as.data.frame(cbind(pct_signal=pct_signal, pv=pv_min, patterns_to_test[,c('r0', 'r1')])), 
       args=list(n=n, m=m, size_pattern=size_pattern, max_size_pattern_to_test=max_size_pattern_to_test))
}
