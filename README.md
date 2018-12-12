# Adaptive Rejection Sampling
R Package: Adaptive Reject Sampling (`ars`)

`ars` can quickly generates number of observations sampling from log-concave density function that is hard to evaluate using adaptive rejection sampling. For more details about the algorithms of adaptive rejection sampling, see Gilks et al. (1992). The implemented method is tangent approach.


## Download the latest version
Development of `ars` can be tracked at https://github.com/yantingpan/ars.
Package can be installed via devtools::install_github('yantingpan/ars')


## Performance of sampling
Several log-concave density functions(trunacted normal, gamma, beta, etc.) are used to test the performance of `ars`. The density of resulting sample points genrated by `ars` fits well with the input density function.   


## Authors
Vincent Myers <vincent_myers@berkeley.edu>, Yanting Pan <yanting_pan@berkeley.edu>, Zhenni Ye <ye.zhenni@berkeley.edu>.
[For any mistakes found, please accept our apology and email Yanting through <yanting_pan@berkeley.edu> so we can fix it as soon as possible.]
