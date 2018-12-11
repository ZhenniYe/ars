# Adaptive Rejection Sampling
R Package: Adaptive Reject Sampling (`ars`)

`ars` can quickly generates number of observations with density function that is hard to evaluate. For more details about algorithms of rejection sampling, see Gilks et al. (1992). The implemented method is tangent approach instead of secant approach.


## Download the latest version
Development of `ars` can be tracked at https://github.com/yantingpan/ars.
Package can be installed by using devtools::install_github(’yantingpan/ars’)


## Performance of sampling
Several log-concave density functions(trunacted normal, gamma, beta, etc.) are used to test the performance of `ars`. The density of sample points genrated by `ars` is nearly the same as the sampled density function.   


## Contributors
This package is written by Yanting Pan, Vince Mayers and Zhenni Ye. If mistake is found, please email Yanting through <yanting_pan@berkeley.edu>, and we will fix the problems as soon as possible.
