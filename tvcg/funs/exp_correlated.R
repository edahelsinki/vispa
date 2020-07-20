source("null_distributions.R")
source("pvalue.R")
library(stats)
library(ggplot2)
library(tidyr)
library(dplyr)

#' Experiment with correlated test statistics
#' 
#' Demonstrate that when test statistics are correlated, minp has more power 
#' than other methods. Create correlated test statistics (using Gaussian processes), 
#' "plant" an outlier point, then compute a p-value for the outlier using various 
#' methods.
#'
#' @param n_grid vector of numbers of test statistics
#' @param l_grid vector of Gaussian process length scales
#' @param outlier_value 
#' @param n_resample number of samples to use for p-value
#' @param pv_funs list of functions f(t0, t_null) that output p=values for t0
#'
#' @return p-value of outlier
#' @export
#'
#' @examples
#' set.seed(2020)
#' n_grid = c(10, 1000)
#' l_grid = c(10, 100, 1000)
#' exp1 = do_exp_correlated(n_grid, l_grid)
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
      cat(paste0("Computing using parameters ", i, "/", nrow(params), "\n"))
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
#' @return
#' @export
#'
#' @examples
make_sample_test_stat <- function(n, l) {
  gp = train_gp(1:n, kern = function(x,y) kernel_sqexp(x,y,l=l))
  function(n=1) gp$prior(n)
}

#' Plot examples of correlated test statistics
#'
#' @param n 
#' @param l 
#' @param outlier_value 
#'
#' @return
#' @export
#'
#' @examples
plot_example = function(n, l, outlier_value) {
  n_null = 10
  sample_test_stat = make_sample_test_stat(n, l)  
  t0 = sample_test_stat()
  t0[1] = outlier_value
  t_null = sample_test_stat(n_null)
  
  df = data.frame(i=1:n, 
                  x=t0)
  df_surrogates = t_null %>% 
    as.data.frame() %>% 
    dplyr::mutate(i=1:n) %>% 
    tidyr::pivot_longer(-i)
  gg_surrogates = geom_line(data=df_surrogates, aes(i, value, group=name), col="grey", alpha=0.5)
  ggplot(df) + 
    gg_surrogates + 
    geom_line(aes(i, x), lwd=0.5) + 
    theme_classic() + 
    theme(axis.line.y = element_blank(), 
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank()) +
    labs(y="", x="") 
}
