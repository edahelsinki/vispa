required_packages = c("dplyr",
                      "ggplot2",
                      "ggpubr",
                      "lubridate",
                      "rmarkdown",
                      "scagnostics",
                      "sp",
                      "stats",
                      "tidyr")

packages_to_install = required_packages[!required_packages %in% installed.packages()]
for (pkg in packages_to_install) install.packages(pkg)
