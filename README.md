
<!-- README.md is generated from README.Rmd. Please edit that file -->
robsurvey
=========

The package **robsurvey** provides several functions to compute robust survey statistics. The package supports the computations of robust means, totals, and ratios. Available methods are Huber M-estimators, trimming, and winsorization. **robsurvey** complements the famous [survey package](https://cran.r-project.org/web/packages/survey/index.html).

Overview
--------

The following functions are provided in **robsurvey**:

-   `weighted.median()`
-   `weighted.quantile()`
-   `weighted.mad()`
-   `weighted.total()`
-   `weighted.mean()`
-   `msvymean()`
-   `msvytotal()`
-   `tsvymean()`
-   `tsvytotal()`
-   `wsvymean()`
-   `wsvytotal()`

Installation
------------

You can install robsurvey from github with:

``` r
# install.packages("devtools")
devtools::install_github("martinSter/robsurvey")
```

Example
-------

In the following example, we showcase a typical use of the package **robsurvey**. The data we use are from the package **survey** and describe the student performance in California schools. We will show different ways of how to compute a robust mean value for the Academic Performance Index (API) in 2000. The variable is denoted as `api00`. The following code chunk simply loads the data and defines the survey design (based on the **survey** package).

``` r
## basic example code
library(robsurvey)
#> Loading required package: survey
#> Warning: package 'survey' was built under R version 3.5.3
#> Loading required package: grid
#> Loading required package: Matrix
#> Loading required package: survival
#> 
#> Attaching package: 'survey'
#> The following object is masked from 'package:graphics':
#> 
#>     dotchart
#> 
#> Attaching package: 'robsurvey'
#> The following object is masked from 'package:stats':
#> 
#>     weighted.mean

# load the api dataset
data(api)

# define survey design
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, 
                    data = apistrat, fpc = ~fpc)
```

In the following code chunk we first compute the robust Horvitz-Thompson M-estimator for the mean. In addition, we compute the trimmed and winsorized (robust) Horvitz-Thompson mean. Note how the estimates and the corresponding standard errors vary.

``` r
# compute the robust Horvitz-Thompson M-estimator of the mean
msvymean(~api00, dstrat, k = 2)
#>          mean    SE
#> api00 662.907 8.926

# compute the robust trimmed Horvitz-Thompson mean
tsvymean(~api00, dstrat, k = 2)
#>          mean     SE
#> api00 655.362 11.568

# compute the robust winsorized Horvitz-Thompson mean
wsvymean(~api00, dstrat, k = 2)
#>          mean     SE
#> api00 640.599 11.568
```

It is also possible to use `msvymean()` in combination with `svyby()` from the **survey** package. The variable `stype` denotes the school level: elementary, middle, and high school.

``` r
# Domain estimates
svyby(~api00, by = ~stype, design = dstrat, msvymean, k = 1.34)
#>   stype    api00       se
#> E     E 675.8203  9.44767
#> H     H 629.1850 10.88119
#> M     M 635.1765 12.63996
```
