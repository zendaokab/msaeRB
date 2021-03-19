
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

## Example

Example for univariate:

``` r
library(msaeRB)
# Using parameter 'data'
data("datamsaeRB")

# Estimation
est_sae = est_saeRB(Y1 ~ X1 + X2, v1, w1, data = datamsaeRB)

# MSE
mse_sae = mse_saeRB(Y1 ~ X1 + X2, v1, w1, data = datamsaeRB)
#> 
#> Bootstrap procedure with B = 1000 iterations starts.
#> 
#> Bootstrap procedure with B = 1000 iterations ends.

# Show the result of est_sae
est_sae
#> $eblup
#> $eblup$est.eblup
#>          Y1
#> 1  11.63062
#> 2  11.85871
#> 3  12.37769
#> 4  12.00198
#> 5  11.48357
#> 6  11.55640
#> 7  12.33156
#> 8  12.44401
#> 9  11.91029
#> 10 12.39311
#> 11 12.19242
#> 12 11.92161
#> 13 11.98016
#> 14 11.52225
#> 15 11.81784
#> 16 12.27088
#> 17 12.44576
#> 18 11.73423
#> 19 12.48320
#> 20 12.33839
#> 21 11.79153
#> 22 12.34885
#> 23 11.66482
#> 24 11.75552
#> 25 11.55350
#> 26 12.26371
#> 27 12.31559
#> 28 12.24639
#> 29 12.23320
#> 30 12.50452
#> 
#> $eblup$est.eblupRB
#>          Y1
#> 1  11.63077
#> 2  11.85885
#> 3  12.37784
#> 4  12.00213
#> 5  11.48371
#> 6  11.55654
#> 7  12.33171
#> 8  12.44416
#> 9  11.91043
#> 10 12.39326
#> 11 12.19257
#> 12 11.92175
#> 13 11.98030
#> 14 11.52239
#> 15 11.81798
#> 16 12.27102
#> 17 12.44591
#> 18 11.73437
#> 19 12.48335
#> 20 12.33854
#> 21 11.79167
#> 22 12.34900
#> 23 11.66496
#> 24 11.75566
#> 25 11.55364
#> 26 12.26386
#> 27 12.31574
#> 28 12.24654
#> 29 12.23335
#> 30 12.50467
#> 
#> 
#> $fit
#> $fit$method
#> [1] "REML"
#> 
#> $fit$convergence
#> [1] TRUE
#> 
#> $fit$iteration
#> [1] 2
#> 
#> $fit$estcoef
#>                   beta std. error   t value      p-value
#> (Intercept)  9.9210029 1.40048406  7.083981 1.400708e-12
#> X1          -0.2513358 0.03487236 -7.207307 5.706928e-13
#> X2           0.4605143 0.12968594  3.550996 3.837763e-04
#> 
#> $fit$refvar
#>            [,1]
#> [1,] 0.02597105
#> 
#> 
#> $agregation
#>                              Y1
#> agregation.direct      12.03692
#> agregation.eblup       12.03677
#> agregation.eblup.ratio 12.03692
```

Example for multivariate:

``` r
library(msaeRB)
# Without parameter 'data'
Fo = list(f1 = datamsaeRB$Y1 ~ datamsaeRB$X1 + datamsaeRB$X2,
          f2 = datamsaeRB$Y2 ~ datamsaeRB$X1 + datamsaeRB$X2,
          f3 = datamsaeRB$Y3 ~ datamsaeRB$X1 + datamsaeRB$X2)
vardir = datamsaeRB[, c("v1", "v12", "v13", "v2", "v23", "v3")]
weight = datamsaeRB[, c("w1", "w2", "w3")]

# Estimation
est_msae = est_msaeRB(Fo, vardir, weight)

# MSE
mse_msae = mse_msaeRB(Fo, vardir, weight)
#> 
#> Bootstrap procedure with B = 1000 iterations starts.
#> 
#> Bootstrap procedure with B = 1000 iterations ends.

# Show the result of MSE
mse_msae
#> $mse.eblup
#>    datamsaeRB$Y1 datamsaeRB$Y2 datamsaeRB$Y3
#> 1    0.007010474    0.01272565    0.01829645
#> 2    0.007397799    0.01367314    0.01982446
#> 3    0.007277975    0.01338002    0.01935175
#> 4    0.007024893    0.01276092    0.01835333
#> 5    0.007525790    0.01398624    0.02032938
#> 6    0.007271614    0.01336446    0.01932665
#> 7    0.007106428    0.01296037    0.01867499
#> 8    0.007227220    0.01325586    0.01915152
#> 9    0.007144129    0.01305260    0.01882372
#> 10   0.007185658    0.01315419    0.01898756
#> 11   0.007150445    0.01306805    0.01884864
#> 12   0.007010223    0.01272503    0.01829546
#> 13   0.007233692    0.01327169    0.01917705
#> 14   0.007271281    0.01336364    0.01932534
#> 15   0.007154048    0.01307686    0.01886285
#> 16   0.007293023    0.01341683    0.01941111
#> 17   0.007225104    0.01325069    0.01914317
#> 18   0.007102008    0.01294956    0.01865755
#> 19   0.007317367    0.01347638    0.01950715
#> 20   0.007883370    0.01486097    0.02174005
#> 21   0.007403115    0.01368614    0.01984543
#> 22   0.007338050    0.01352698    0.01958875
#> 23   0.007245573    0.01330076    0.01922392
#> 24   0.007313363    0.01346659    0.01949136
#> 25   0.007443766    0.01378559    0.02000580
#> 26   0.007116563    0.01298517    0.01871497
#> 27   0.007074608    0.01288253    0.01854946
#> 28   0.007025752    0.01276302    0.01835672
#> 29   0.007076076    0.01288612    0.01855525
#> 30   0.007288537    0.01340586    0.01939341
#> 
#> $pbmse.eblupRB
#>    datamsaeRB$Y1 datamsaeRB$Y2 datamsaeRB$Y3
#> 1    0.007015449    0.01280613    0.01848852
#> 2    0.007377658    0.01366118    0.01990889
#> 3    0.007288419    0.01338543    0.01944527
#> 4    0.007029372    0.01276679    0.01851608
#> 5    0.007499085    0.01387073    0.02025959
#> 6    0.007284548    0.01345583    0.01940157
#> 7    0.007112686    0.01298647    0.01886326
#> 8    0.007227051    0.01328476    0.01927360
#> 9    0.007151277    0.01312528    0.01902267
#> 10   0.007173481    0.01321076    0.01905874
#> 11   0.007157801    0.01309861    0.01901040
#> 12   0.007018966    0.01276308    0.01857168
#> 13   0.007238103    0.01331737    0.01933909
#> 14   0.007250264    0.01338316    0.01935921
#> 15   0.007147256    0.01305434    0.01894990
#> 16   0.007283234    0.01339014    0.01956173
#> 17   0.007241212    0.01331154    0.01927470
#> 18   0.007118923    0.01299011    0.01878413
#> 19   0.007289666    0.01337688    0.01947539
#> 20   0.007800953    0.01465303    0.02146016
#> 21   0.007381786    0.01362265    0.01988992
#> 22   0.007321924    0.01348304    0.01969166
#> 23   0.007244528    0.01328253    0.01917915
#> 24   0.007317799    0.01345850    0.01956893
#> 25   0.007412795    0.01370840    0.01996931
#> 26   0.007127487    0.01301122    0.01883017
#> 27   0.007079342    0.01294654    0.01869917
#> 28   0.007045588    0.01283910    0.01859037
#> 29   0.007074784    0.01287507    0.01865532
#> 30   0.007254225    0.01337512    0.01953320
#> 
#> $running.time
#> Time difference of 4.60364 mins
```
