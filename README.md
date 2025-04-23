# Resampling method for density estimation at a survival quantile
We propose a resampling procedure that allows to estimate $f(F^{-1}(p))$ the survival density at a given quantile, for a given probability $0 < p < 1$.
We present a method inspired by Lin et al., 2015 in order to estimate the density $f$ evaluated at the quantile $F^{-1}(p)$.
We require the densities at the quantiles to be strictly positive, and we denote as $\hat{F}$ the consistent estimator for $F$ obtained from the usual Kaplan-Meier estimation. Taking this estimator we obtain $\hat{F}^{-1}(p)$ the estimators of the inverse distribution at $p$. Then we propose the following estimation procedure:

1. Generate B realizations of the gaussian $T \sim \mathcal{N}(0, \sigma^2)$, denoted by $T_1,..., T_B$
2.  Calculate $\sqrt{n}\left( \hat{F}\left(\hat{F}^{-1}(p) + \dfrac{T_b}{\sqrt{n}}\right) - p \right), b = 1,..., B$ and denote them as $y_b$, then the least squares estimate of $f(F^{-1}(p))$ is $\hat{A} = (x'x)^{-1}x'Y$, where $x= (T_1,..., T_B)^T$ and $Y = (y_1,..., y_B)^T$

In order to choose the variance of the generated gaussians $\sigma^2$, we suggest to perform an automatic bandwidth selection via grid-search.

## Variance selection 
In our procedure, we need to specify the variance of the generated gaussian variables. Analysing the MSE for several distributions (such as exponential, piecewise exponential, cauchy, weibull...), we notice that, for small sample sizes, the value of $\sigma^2$ plays an important role in the value of the MSE. Indeed, as one increases sample size, the region where MSE is minimized becomes larger, which allows us to choose the variance more freely in a large interval of values where the MSE is small.

In the simulations performed to study the behavior of the MSE for different values of $\sigma^2$, we notice for all sample sizes that the MSE decreases until it reaches a plateau where it stays small for an interval of values of $\sigma$, after which it increases again. Our goal is then to select a value of $\sigma$ that is inside such interval, which grows larger as sample size increases. 

The behavior of the MSE for different values of $\sigma$ is illustrated below, for 100 repetitions of the code for an exponential distribution with rate 1.5 evaluated at the median, with exponential censoring with rate 0.48 and B = 100000 generated gaussians.

![new_expo_n50_M100](https://github.com/user-attachments/assets/572cc906-e562-4f06-8847-dc4873499e58)

![new_expo_n200_M100](https://github.com/user-attachments/assets/7c72629f-3ff9-434e-a1b8-ca12dc407801)

![new_expo_n1000_M100](https://github.com/user-attachments/assets/ca0f7a68-b66a-423b-9898-464795b7f8b8)

We propose a grid-search algorithm in order to choose a value of $\sigma$ in an automatic way. Seen that the best values of this parameter knowing the true value of the density at the quantile belong to an interval, we aim to detect this plateau from the estimated densities at the quantiles, and choose a value that is inside this interval.

We illustrate this procedure for an exponential distribution with rate 1.5 with censoring that follows an exponential with rate 0.12 (which gives approximately 10% of censored observations). We are interested in estimating the density of the exponential at the median, which has the real value equal to 0.75. We then perform the estimation of the density at the quantile for a grid of $\sigma$, and we observe that there exists an interval of values of $\sigma$ such that the estimation is small before it and decreases after it as well. We aim to choose a value of $\sigma$ that lays inside this interval, because this is the interval that corresponds to the values close to the real value (marked in the dotted red line). Our code allows us to detect such interval and then it chooses a value of $\sigma$ that lays inside it. The chosen value is marked in blue in the corresponding figures. In this illustrative example, we obtain an estimation of the density at the median equal to 0.745831 (where the real value is 0.75), for a corresponding value of $\sigma$ equal to 3.75. This is performed for a grid of $\sigma$ that goes from 0.1 to 10, in steps of length 0.05. 

![example of density estim at exponential](https://github.com/user-attachments/assets/51a4f61a-6b2e-4cbf-8f38-0a4ad79a8388)

![diff for density estim exponential](https://github.com/user-attachments/assets/2c19be49-ae6a-4f97-8454-845116ead13f)

All simulations can be run using the auximilary codes in src/code_source and src/simulations_and_results
The folder src has the corresponding codes for computing the estimation of the density for the Exponential and Cauchy distributions. The estimations using our proposed procedure are compared to the ones obtained using kernel density estimation with bandwidth parameter approximated using leave-one-out cross-validation.

This repository is organised as follows:
- /src:
    - code_source: auxiliary codes for test in presence of data and for planning clinical trial using the kernel density estimation and our resmpling procedure
    - planning a clinical trial: codes for the section of Planning a clinical trial from our paper "Univariate and multivariate test of equality of quantiles with right-censored data"
    - simulations_and_results: comparing density estimation with kernel density and our resampling method
    - test on OAK data: applying the density estimatoin for the univariate and multivariate test of equality of quantiles detailed in our paper "Univariate and multivariate test of equality of quantiles with right-censored data"
 
- /MSE: codes for the figures for MSE varying with variance, for exponential distribution

# References
Lin, C., Zhang, L., & Zhou, Y. (2015). Conditional quantile residual lifetime models for right censored data. Lifetime data analysis, 21, 75-96.
