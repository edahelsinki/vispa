library(stats)
library(tidyr)
library(dplyr)
source("funs/null_distributions.R")
source("funs/pvalue.R")

#' Experiment with correlated test statistics
#' 
#' Create correlated test statistics (using Gaussian processes), "plant" an 
#' outlier point, then compute a p-value for the outlier using various methods.
#'
#' @param n_grid vector of numbers of test statistics
#' @param l_grid vector of Gaussian process length scales
#' @param outlier_value 
#' @param n_resample number of samples to use for p-value
#' @param pv_funs list of functions f(t0, t_null) that output p=values for t0
#'
#' @return p-value of outlier
do_exp_correlated <- function(n_grid = c(10, 1000), 
                              l_grid = c(3, 10, 100), 
                              outlier_value = 3, 
                              n_resample = 1000, 
                              n_rep = 10,
                              pv_funs = list(raw=pvalue_vec, 
                                             minp=minP, 
                                             bonf=function(t0, t_null) p.adjust(pvalue_vec(t0, t_null)), 
                                             BH=function(t0, t_null) p.adjust(pvalue_vec(t0, t_null), "BH"), 
                                             BY=function(t0, t_null) p.adjust(pvalue_vec(t0, t_null), "BY")
                              ), 
                              verbose=FALSE) {
  params = expand.grid(n = n_grid, l = l_grid, KEEP.OUT.ATTRS = FALSE) 
  pv = array(0, c(nrow(params), length(pv_funs), n_rep), list(NULL, names(pv_funs), NULL))
  for (i in 1:nrow(params)) {
    if (verbose)
      cat(paste0("Parameter set ", i, "/", nrow(params), "\n"))
    sample_test_stat = make_sample_test_stat(params$n[i], params$l[i])  
    for (k in 1:n_rep) {
      if (verbose && (k %% 10 == 0)|k==1)
        cat(paste0("\trep", k, "/", n_rep, "\n"))
      t0 = sample_test_stat()
      t0[1] = outlier_value
      t_null = sample_test_stat(n_resample)  
      for (j in seq_along(pv_funs)) {
        # Save only the p-value of the outlier
        pv[i, j, k] = pv_funs[[j]](t0, t_null)[1]
      }
    }
  }
  pv_avg = apply(pv, c(1,2), median)
  list(result = cbind(params, pv = pv_avg),
       pv_array = pv,
       args = list(n_grid = n_grid,
                   l_grid = l_grid,
                   outlier_value = outlier_value,
                   n_resample = n_resample, 
                   n_rep = n_rep,
                   pv_funs = pv_funs))
  cbind(params, pv = pv_avg)
}

#' Wrapper for making a Gaussian process prior
#'
#' @param n 
#' @param l 
#'
#' @return function
make_sample_test_stat <- function(n, l) {
  gp = train_gp(1:n, kern = function(x,y) kernel_sqexp(x,y,l=l))
  function(n=1) gp$prior(n)
}


# Run simulation ####

filename = "data/exp_correlated.rds"
if (file.exists(filename)) {
  exp_correlated = readRDS(filename)
} else {
  set.seed(2020)
  n_rep = 100
  exp_correlated = do_exp_correlated(n_grid = c(10, 1000), 
                                     l_grid = c(3, 10, 25, 100, 250, 500, 1000), 
                                     n_rep = n_rep, 
                                     verbose = TRUE)
  saveRDS(exp_correlated, filename)
}
