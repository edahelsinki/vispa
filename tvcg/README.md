# README

To reproduce all results run: 
```bash
Rscript --vanilla make.R
```

Note that a ~200Mb dataset will be downloaded into `data/banana-data.rds`.
The experiments can also be run individually as `Rscript --vanilla run_exp_*.R`. 
The results are stored in `results/`. To reproduce the results yourself, delete the files in `results/`. 

Alternatively, open RStudio and run `tvcg_examples.Rmd` and `tvcg_experiments.Rmd`.
