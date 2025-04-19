
<!-- README.md is generated from README.Rmd. Please edit that file -->

# Max-Type Test for Marginal Correlation with Bootstrap (MarginalMaxTest)

<!-- badges: start -->
<!-- badges: end -->

The goal of MarginalMaxTest is to test the marginal correlation between
a scalar response variable with a vector of explanatory variables using
the max-type test with bootstrap. The test is based on the max-type
statistic and its asymptotic distribution under the null hypothesis of
no marginal correlation. The bootstrap procedure is used to approximate
the null distribution of the test statistic. The package provides a
function for performing the test. For more technical details, refer to
Zhang and Laber (2014) <doi:10.1080/01621459.2015.1106403>.

## Installation

Install the R CRAN version of MarginalMaxTest like so:

``` r
install.packages("MarginalMaxTest")
```

You can also install the development version of MarginalMaxTest from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("canyi-chen/MarginalMaxTest")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(MarginalMaxTest)
# Generate sample data
set.seed(47)
n <- 200
p <- 10
x <- matrix(rnorm(n*p), n, p)
y <- 0.25*x[,1] + rnorm(n)
# Run the test
marginal.test(x, y, B = 200, method = "adaptive")
#> $p_value
#> [1] 0.01
#> 
#> $time
#> [1] 0.007
marginal.test(x, y, B = 200, method = "max")
#> $p_value
#> [1] 0
#> 
#> $time
#> [1] 0.007
marginal.test(x, y, B = 200, method = "sum")
#> $p_value
#> [1] 0.11
#> 
#> $time
#> [1] 0.007
```
