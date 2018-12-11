# Adaptive Rejection Sampling
R Package: Adaptive Reject Sampling (`ars`)

`ars` can quickly generates number of observations sampling from log-concave density function that is hard to evaluate using adaptive rejection sampling. For more details about the algorithms of adaptive rejection sampling, see Gilks et al. (1992). The implemented method is tangent approach.


## Download the latest version
Development of `ars` can be tracked at https://github.com/yantingpan/ars.
Package can be installed via devtools::install_github('yantingpan/ars')


## Performance of sampling
Several log-concave density functions(trunacted normal, gamma, beta, etc.) are used to test the performance of `ars`. The density of resulting sample points genrated by `ars` fit well with the input density function.   


## Authors
Myers, Vince <vincent_myers@berkeley.edu>; Pan, Yanting <yanting_pan@berkeley.edu>; Ye, Zhenni <ye.zhenni@berkeley.edu>;  based on Gilks and Wild (1992).
[For any mistakes found, please accept our apology and email Yanting through <yanting_pan@berkeley.edu> so we can fix it as soon as possible.]
