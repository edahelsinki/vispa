required_packages = c("dplyr",
                      "forcats",
                      "ggplot2",
                      "ggpubr",
                      "lubridate",
                      "onlineFDR",
                      "purrr",
                      "readr",
                      "scagnostics",
                      "sp",
                      "stats",
                      "stringr",
                      "tidyr")

packages_to_install = required_packages[!required_packages %in% installed.packages()]
for (pkg in packages_to_install) install.packages(pkg)
