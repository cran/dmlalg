
<!-- README.md is generated from README.Rmd. Please edit that file -->

# dmlalg

The *dmlalg* package contains implementations of double machine learning
algorithms in *R*. It includes the following functionality:

## Partially linear models with confounding variables

This first set of functions perform inference for the linear parameter
in partially linear models with confounding variables. The standard DML
estimator of the linear parameter has a two-stage least squares
interpretation, which can lead to a large variance and overwide
confidence intervals. We apply regularization to reduce the variance of
the estimator, which produces narrower confidence intervals that remain
approximately valid. Nuisance terms can be flexibly estimated with
machine learning algorithms.

This method is described in Emmenegger and Bühlmann (2021) and
implemented in the function *regsdml*.

  - *regsdml* computes the estimate of the linear parameter in a
    partially linear model with endogenous variables with regularized
    and standard double machine learning methods.
  - *summary.regsdml* is a *summary* method for objects fitted with
    *regsdml*.
  - *confint.regsdml* is a *confint* method for objects fitted with
    *regsdml*.
  - *coef.regsdml* is a *coef* method for objects fitted with *regsdml*.
  - *vcov.regsdml* is a *vcov* method for objects fitted with *regsdml*.
  - *print.regsdml* is a *print* method for objects fitted with
    *regsdml*.

## Installation

You can install the released version of dmlalg from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("dmlalg")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(dmlalg)

## Generate some data:
set.seed(19)
# true linear parameter
beta0 <- 1
n <- 40
# observed confounder
w <- pi * runif(n, -1, 1)
# instrument
a <- 3 * tanh(2 * w) + rnorm(n, 0, 1)
# unobserved confounder
h <- 2 * sin(w) + rnorm(n, 0, 1)
# linear covariate
x <- -1 * abs(a) - h - 2 * tanh(w) + rnorm(n, 0, 1)
# response
y <- beta0 * x - 3 * cos(pi * 0.25 * h) + 0.5 * w ^ 2 + rnorm(n, 0, 1)

## Estimate the linear coefficient from x to y
## (The parameters are chosen small enough to make estimation fast):
## Caveat: A spline estimator is extrapolated, which raises a warning message. 
## Extrapolation lies in the nature of our method. To omit the warning message
## resulting from the spline estimator, another estimator may be used. 
fit <- regsdml(a, w, x, y,
               gamma = exp(seq(-4, 1, length.out = 4)),
               S = 3,
               do_regDML_all_gamma = TRUE,
               cond_method = c("forest",  # for E[A|W]
                               "spline",  # for E[X|W]
                               "spline"), # for E[Y|W]
               params = list(list(ntree = 1), NULL, NULL))
#> Warning in regsdml(a, w, x, y, gamma = exp(seq(-4, 1, length.out = 4)), : 
#> Warning messages:
#> some 'x' values beyond boundary knots may cause ill-conditioned bases
## parm = c(2, 3) prints an additional summary for the 2nd and 3rd gamma-values
summary(fit, parm = c(2, 3),
        correlation = TRUE,
        print_gamma = TRUE) 
#> 
#> Coefficients :
#> regsDML (2.72e+00) :
#>    Estimate Std. Error  z value     Pr(>|z|)
#> b1 0.910255  0.1731559 5.256852 1.465421e-07
#> 
#> regDMLall (9.70e-02) :
#>     Estimate Std. Error  z value     Pr(>|z|)
#> b1 0.7986392  0.1514027 5.274935 1.328031e-07
#> 
#> regDMLall (5.13e-01) :
#>    Estimate Std. Error  z value     Pr(>|z|)
#> b1 0.846176  0.1651298 5.124308 2.986318e-07
#> 
#> 
#> Variance-covariance matrices :
#> regsDML (2.72e+00) :
#>            b1
#> b1 0.02998297
#> 
#> regDMLall (9.70e-02) :
#>            b1
#> b1 0.02292277
#> 
#> regDMLall (5.13e-01) :
#>            b1
#> b1 0.02726785
confint(fit, parm = c(2, 3),
        print_gamma = TRUE) 
#> Two-sided confidence intervals at level 0.95 : 
#> 
#> regsDML (2.72e+00) :
#>        2.5 %   97.5 %
#> b1 0.5708757 1.249634
#> 
#> regDMLall (9.70e-02) :
#>        2.5 %   97.5 %
#> b1 0.5018955 1.095383
#> 
#> regDMLall (5.13e-01) :
#>        2.5 %   97.5 %
#> b1 0.5225276 1.169824
coef(fit) # coefficients
#>     regsDML
#> b1 0.910255
vcov(fit) # variance-covariance matrices
#> 
#> Variance-covariance matrices :
#> regsDML :
#>            b1
#> b1 0.02998297

## Alternatively, provide the data in a single data frame
## (see also caveat above):
data <- data.frame(a = a, w = w, x = x, y = y)
fit <- regsdml(a = "a", w = "w", x = "x", y = "y", data = data,
               gamma = exp(seq(-4, 1, length.out = 4)),
               S = 3)
#> Warning in regsdml(a = "a", w = "w", x = "x", y = "y", data = data, gamma = exp(seq(-4, : 
#> Warning messages:
#> some 'x' values beyond boundary knots may cause ill-conditioned bases

## With more realistic parameter choices:
if (FALSE) {
  fit <- regsdml(a, w, x, y,
                 cond_method = c("forest",  # for E[A|W]
                                 "spline",  # for E[X|W]
                                 "spline")) # for E[Y|W]
  summary(fit)
  confint(fit)

  ## Alternatively, provide the data in a single data frame:
  ## (see also caveat above):
  data <- data.frame(a = a, w = w, x = x, y = y)
  fit <- regsdml(a = "a", w = "w", x = "x", y = "y", data = data)
}
```

<div id="refs" class="references">

<div id="ref-Emmenegger2021">

Emmenegger, Corinne, and Peter Bühlmann. 2021. “Regularizing Double
Machine Learning in Partially Linear Endogenous Models.” *Preprint
arXiv:2101.12525*. <https://arxiv.org/abs/2101.12525>.

</div>

</div>
