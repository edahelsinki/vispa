# Install dependencies
cat("Checking dependencies...\n")
source("install_dependencies.R")

# Download data
cat("Checking SMEAR data...\n")
source("funs/read_data.R")
check_smear()

# Knit .Rmd files
cat("Creating reports...\n")
suppressPackageStartupMessages({
  rmarkdown::render("use_cases.Rmd")
  rmarkdown::render("schematic_examples.Rmd")
})
cat("Done.")
