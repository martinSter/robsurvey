
<!-- README.md is generated from README.Rmd. Please edit that file -->
robsurvey
=========

The package **robsurvey** provides several functions to compute robust survey statistics. The package supports the computations of robust means, totals, and ratios. Available methods are Huber M-estimators, trimming, and winsorization. **robsurvey** complements the famous [survey package](https://cran.r-project.org/web/packages/survey/index.html).

Overview
--------

The following functions are provided in **robsurvey**:

-   `weighted_median()`
-   `weighted_quantile()`
-   `weighted_mad()`
-   `weighted_total()`
-   `weighted_mean()`
-   `svymean_huber()`
-   `svytotal_huber()`
-   `svymean_trimmed()`
-   `svytotal_trimmed()`
-   `svymean_winsorized()`
-   `svytotal_winsorized()`
-   `weighted_line()`
-   `weighted_median_line()`
-   `weighted_median_ratio()`

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
# load and attach packages
library(robsurvey)
library(survey)

# load the api dataset
data(api)

# define survey design
dstrat <- svydesign(id = ~1, strata = ~stype, weights = ~pw, 
                    data = apistrat, fpc = ~fpc)
```

In the following code chunk we first compute the robust Horvitz-Thompson M-estimator for the mean. In addition, we compute the trimmed and winsorized (robust) Horvitz-Thompson mean. Note how the estimates and the corresponding standard errors vary. The scale estimate used in the Huber-type robust M-estimator is the MAD, which is rescaled to be consistent at the normal distribution, i.e. multiplied by the constant 1.4826. The default tuning constant of the Huber-type Horvitz-Thompson M-estimator is *k* = 1.5.

``` r
# compute the robust Horvitz-Thompson M-estimator of the mean
svymean_huber(~api00, dstrat, k = 2)
#>          mean    SE
#> api00 662.907 8.926

# compute the robust trimmed Horvitz-Thompson mean
svymean_trimmed(~api00, dstrat, k = 2)
#>          mean     SE
#> api00 655.362 11.568

# compute the robust winsorized Horvitz-Thompson mean
svymean_winsorized(~api00, dstrat, k = 2)
#>          mean     SE
#> api00 640.599 11.568
```

It is also possible to use `svymean_huber()` in combination with `svyby()` from the **survey** package. The variable `stype` denotes the school level: elementary, middle, and high school.

``` r
# Domain estimates
svyby(~api00, by = ~stype, design = dstrat, svymean_huber, k = 1.34)
#>   stype    api00       se
#> E     E 675.8203  9.44767
#> H     H 629.1850 10.88119
#> M     M 635.1765 12.63996
```

For simulations and as intermediate results the above functions can also be used without the **survey** package. They deliver only the bare estimate.

``` r
# bare-bone function to compute robust Horvitz-Thompson M-estimator of the mean
weighted_mean_huber(apistrat$api00, weights(dstrat), k = 2)
#> [1] 662.9068
```
