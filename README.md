
[![Travis build
status](https://travis-ci.com/zendaokab/msaeRB.svg?branch=main)](https://travis-ci.com/zendaokab/msaeRB)
[![AppVeyor build
status](https://ci.appveyor.com/api/projects/status/github/zendaokab/msaeRB?branch=main&svg=true)](https://ci.appveyor.com/project/zendaokab/msaeRB)
[![Codecov test
coverage](https://codecov.io/gh/zendaokab/msaeRB/branch/main/graph/badge.svg)](https://codecov.io/gh/zendaokab/msaeRB?branch=main)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# msaeRB

Implements multivariate ratio benchmarking Small Area Estimation. This
package provides ratio benchmarking estimation for univariate and
multivariate Small Area Estimation and its MSE. In fact, MSE estimators
for ratio renchmark are not readily available, so resampling method that
called parametric bootstrap is applied. The ratio benchmark model and
parametric bootstrap in this package are based on the model proposed in
Small Area Estimaton (J.N.K Rao and Isabel Molina, 2015)

## Installation

You can install the released version of msaeRB from
[CRAN](https://CRAN.R-project.org) with:

``` r
install.packages("msaeRB")
```

## Atuhors

Zenda Oka Briantiko, Azka Ubaidillah

## Maintainer

Zenda Oka Briantiko <221710087@stis.ac.id>

## Functions

  - est\_saeRB : Produces EBLUPs Ratio Benchmarking based on a
    Univariate Fay-Herriot (Model 1)
  - mse\_saeRB : Parametric Bootstrap Mean Squared Error Estimators of
    Ratio Benchmarking for Univariate Small Area Estimation
  - est\_msaeRB : Produces EBLUPs Ratio Benchmarking based on a
    Multivariate Fay-Herriot (Model 1)
  - mse\_msaeRB : Parametric Bootstrap Mean Squared Error Estimators of
    Ratio Benchmarking for Multivariate Small Area Estimation

## References

  - Rao, J.N.K & Molina. (2015). Small Area Estimation 2nd Edition. New
    York: John Wiley and Sons, Inc.
  - Benavent, Roberto & Morales, Domingo. (2015). Multivariate
    Fay-Herriot models for small area estimation. Computational
    Statistics and Data Analysis 94 2016 372-390. DOI:
    10.1016/j.csda.2015.07.013.
  - Ubaidillah, Azka et al. (2019). Multivariate Fay-Herriot models for
    small area estimation with application to household consumption per
    capita expenditure in Indonesia. Journal of Applied Statistics.
    46:15. 2845-2861. DOI: 10.1080/02664763.2019.1615420.
  - Datta, G.S., Day, B., and Maiti, T. (1998), Multivariate Bayesian
    Small Area Estimation: An Application to Survey and Satellite Data,
    Sankhy¯a, Series A, 60, 1–19.
  - Wang, J., Fuller, W.A., and Qu, Y. (2008). Small Area Estimation
    Under Restriction, Survey Methodology, 34, 29–36
  - Małgorzata Karolina Krzciu (2017). On the Simulation Study of
    Jackknife and Bootstrap MSE Estimators of a Domain Mean Predictor
    for Fay‑Herriot Model. DOI:
    <http://dx.doi.org/10.18778/0208‑6018.331.11>
