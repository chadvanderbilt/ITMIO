packages <- c("tidyverse",
              "gtsummary",
              "varhandle",
              "readxl",
              "ggpubr",
              "janitor",
              "BiocManager",
              "forestmodel",
              "survminer", 
              "survival",
              "scales",
              "lubridate",
              "stringr",
              "webshot2",
              "data.table",
              "viridis", 
              "devtools")
install.packages(packages)
if(getRversion() >= "4.0.0"){
  if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager", repos = "http://cran.us.r-project.org")
  BiocManager::install(packages)
}
install.packages(packages)

# Install Bioconductor packages if needed.
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("EnhancedVolcano")
BiocManager::install("ComplexHeatmap")
install vdbr from source if needed
devtools::install_git("https://github.mskcc.org/vdblabinternal/vdbr.git")
