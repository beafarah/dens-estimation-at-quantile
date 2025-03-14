# Resampling method for density estimation at a survival quantile
We propose a resampling procedure that allows to estimate $f(F^{-1}(p)$ the survival density at a given quantile, for a given probability $0 < p < 1$.
We present a method inspired by Lin et al., 2015 in order to estimate the density $f$ evaluated at the quantile $F^{-1}(p)$.
We require the densities at the quantiles to be strictly positive, and we denote as $\hat{F}$ the consistent estimator for $F$ obtained from the usual Kaplan-Meier estimation. Taking this estimator we obtain $\hat{F}^{-1}(p)$ the estimators of the inverse distribution at $p$. Then we propose the following estimation procedure:

1. Generate B realizations of the gaussian $T sim \mathcal{N}(0, \sigma^2)$, denoted by $T_1,..., T_B$
2.  Calculate $\sqrt{n}\left( \hat{F}\left(\hat{F}^{-1}(p) + \dfrac{\Tau_b}{\sqrt{n}}\right) - p \right), b = 1,..., B$ and denote them as $y_b$, then the least squares estimate of $A(F^{-1}(p))$ is $\hat{A} = (x'x)^{-1}x'Y$, where $x= (\Tau_1,..., \Tau_B)^T$ and $Y = (y_1,..., y_B)^T$

In order to choose the variance of the generated gaussians $\sigma^2$, we suggest to perform an automatic bandwidth selection via grid-search.

## Variance selection 
In our procedure, we need to specify the variance of the generated gaussian variables. Analysing the MSE for several distributions (such as exponential, piecewise exponential, cauchy, weibull...), we notice that, for small sample sizes, the value of $\sigma^2$ plays an important role in the value of the MSE. Indeed, as one increases sample size, the region where MSE is minimized becomes larger, which allows us to choose the variance more freely in a large interval of values where the MSE is small.

In the simulations performed to study the behavior of the MSE for different values of $\sigma^2$, we notice for all sample sizes that the MSE decreases until it reaches a plateau where it stays small for an interval of values of $\sigma$, after which it increases again. Our goal is then to select a value of $\sigma$ that is inside such interval, which grows larger as sample size increases. 
