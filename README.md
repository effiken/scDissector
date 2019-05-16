# scDissector
scDissector is an exploratory data analysis tool for single-cell RNA-seq data implemented as R shiny app.

Please contact Ephraim.Kenigsberg at mssm edu for additional information.

## Installation

1.  Install R (in case you don't have it already installed):

    Download and install [R](https://cran.r-project.org/).
    
    You may want to install [RStudio](https://www.rstudio.com/) as well.
    
2.  Install devtools package in R (in case you don't have it already installed):

    **install.packages("devtools")**
    
    **library(devtools)**
    
3.  Install scDissector:

    **install_github("effiken/scDissector")**
    
    *or alternatively if you would like to take a risk...*
    
    **install_github("effiken/scDissector",ref = "devel")**
    
    *see the [readme of the devel version](https://github.com/effiken/scDissector/blob/devel/README.md)*

## Update

1. Load devtools

**library(devtools)**

2. Install the package as in (3) above

## Running scDissector in R

**library(scDissector)**

**run_scDissector()**

