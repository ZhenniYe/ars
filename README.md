# Adaptive Rejection Sampling
R Package: Adaptive Reject Sampling (`ars`)

`ars` can quickly generates number of observations with density function that is hard to evaluated. For more details of algorithms with regard to rejection sampling, see Gilks et al. (1992). The method we implemented is tangent approach instead of secant approach. The package is called 'ars' and you can install and use in R. The implemented method is tangent approach instead of secant approach. Package can be installed by using devtools::install_github(’yantingpan/ars’) 


## Download the latest version
Development of `ars` can be tracked at https://github.com/yantingpan/ars .
Package can be installed by using devtools::install_github(’yantingpan/ars’)


## Performance of sampling
Several log-concave density function(trunacted normal, gamma, beta, etc.) are used to test the performance of `ars`. The density of sample points genrated by `ars` is nearly the same as the sampled density function.   


## Contributors
This package is written by Yanting Pan, Vince Mayers, Zhenni Ye. If mistake is found, please emali Yanting through <yanting_pan@berkeley.edu>, and we will fix the problems as soon as possible.
