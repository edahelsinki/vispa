source("funs/exp_overlap.R", chdir = TRUE)
filename = "results/exp_overlap.rds"
overwrite=FALSE
if (file.exists(filename) & overwrite==FALSE) { 
  exp_overlap = readRDS(filename)
} else {
  set.seed(2020)
  exp_overlap = do_experiment_overlap(1000, 100, 100, TRUE)
  saveRDS(exp_overlap, filename)
}
