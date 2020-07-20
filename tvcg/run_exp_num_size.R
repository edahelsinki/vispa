source("funs/exp_num_size.R", chdir=TRUE)
filename = "results/exp_num_size.rds"
overwrite=FALSE
if (file.exists(filename) & overwrite==FALSE) {
  exp_num_size = readRDS(filename)
} else {
  set.seed(2020)
  exp_num_size = do_experiment_num_size(n_grid = 5000, 
                                        m_grid = 5, 
                                        num_tile_grid = c(1, 25, 50, 75, 100), 
                                        size_tile_grid = c(5, 10, 15, 20, 25, 30, 35, 40, 45, 50), 
                                        n_resample = 1000, 
                                        verbose = TRUE)
  saveRDS(exp_num_size, filename)
}
