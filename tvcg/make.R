# Install dependencies
cat("Checking dependencies...\n")
source("install_dependencies.R")

# Download data
cat("Checking smear data...\n")
source("funs/read_data.R")
if (!file.exists("data/banana-data.rds"))
  d=read_smear_all("data/banana-data.rds", "data/DMPS_Event_Classification.txt")

# Run experiments
cat("Running experiments...\n")
suppressPackageStartupMessages({
  source("run_exp_correlated.R", chdir = TRUE)
  source("run_exp_expert_naive.R", chdir = TRUE)
  source("run_exp_num_size.R", chdir = TRUE)
  source("run_exp_overlap.R", chdir = TRUE)
})
# Knit .Rmd files
cat("Creating reports from tvcg_examples.Rmd and tvcg_experiments.Rmd...\n")
suppressPackageStartupMessages({
  rmarkdown::render("tvcg_examples.Rmd")
  rmarkdown::render("tvcg_experiments.Rmd")
})
cat("Done.")
