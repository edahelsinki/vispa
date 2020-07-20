source("funs/exp_correlated.R", chdir = TRUE)

filename = "results/exp_correlated.rds"
if (file.exists(filename)) {
  exp_correlated = readRDS(filename)
} else {
  make_sample_test_stat <- function(n, l) {
    gp = train_gp(1:n, kern = function(x,y) kernel_sqexp(x,y,l=l))
    function(n=1) gp$prior(n)
  }
  set.seed(2020)
  n_rep = 100
  exp_correlated = do_exp_correlated(n_grid = c(10, 1000), 
                                     l_grid = c(3, 10, 25, 100, 250, 500, 1000), 
                                     n_rep = n_rep, 
                                     verbose = TRUE)
  saveRDS(exp_correlated, filename)
}
