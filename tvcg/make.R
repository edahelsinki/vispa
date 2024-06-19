# Install dependencies
cat("Checking dependencies...\n")
source("install_dependencies.R")

# Download data
cat("Checking SMEAR data...\n")
source("funs/read_data.R")
check_smear()

# Run appendix experiment
cat("Running simulation for appendix (~20min)...\n")
source("run_exp_correlated.R")

# Knit .Rmd files
cat("Creating reports...\n")
suppressPackageStartupMessages({
  rmarkdown::render("use_cases.Rmd")
  rmarkdown::render("schematic_examples.Rmd")
  rmarkdown::render("appendix_figures.Rmd")
})
cat("Done.")
