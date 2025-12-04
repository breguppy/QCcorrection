# QCcorrection
Shiny app for QC correction of metabolomics data

### Requirements

 - R >= 4.4.1
 
### Installation

#### Install Via GitHub
To install via GitHub, use the the devtools package. To install and use devtools:
```r
install.packages("devtools")
library(devtools)
```
To install the QCcorrection package:

```r
devtools::install_github("breguppy/QCcorrection")
```

To be able to use all features of the app, please install the following packages:
```r
install.packages(c(
  "randomForest", "robustbase", "outliers", "EnvStats", "ggtext",
  "cowplot", "pagedown", "rmarkdown", "zip", "corpcor",
  "httpuv", "jsonlite", "shinytest2", "testthat", "chromote",
  "pkgload", "knitr", "pagedown", "ggtext"
))
```

### To run App
```r
QCcorrection::run_app()
```

### Example Raw Data Stucture
<img align="center" src="https://github.com/breguppy/QCcorrection/blob/main/www/example_data_structure.png">

# Bug Reports/New Features

#### If you run into any issues or bugs please submit a [GitHub issue](https://github.com/breguppy/QCcorrection/issues) with details of the issue.

- If possible please include a [reproducible example](https://reprex.tidyverse.org/). 

#### Any requests for new features or enhancements can also be submitted as [GitHub issues](https://github.com/breguppy/QCcorrection/issues).

#### [Pull Requests](https://github.com/breguppy/QCcorrection/pulls) are welcome for bug fixes, new features, or enhancements.

## Citation

Coming Soon
