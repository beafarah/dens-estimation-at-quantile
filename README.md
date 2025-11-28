# Resampling method for density estimation at quantile and Clinical trial design for the comparison of single and multiple quantiles with right-censored data
This repository contains functions and code used in the following papers:

1."Designing clinical trials for the comparison of single and multiple quantiles with right-censored data"
2. "A note on a resampling procedure for density estimation at quantiles"

The resampling procedure implemented here was inspired by Lin et al. (2015) and is used to estimate the density at a given quantile $f(F^{-1}(p))$, for a given probability $0 < p < 1$.
This estimation is essential for performing tests of equality of quantiles in the presence of right-censored data, as proposed in our clinical trial design paper.

## Resampling procedure
Let $ \hat{F}$ denote the Kaplan-Meier estimator of the distribution $F$, and $\hat{F}^{-1}(p)$ its inverse at probability $p$. The procedure to estimate $f(\hat{F}^{-1}(p))$ is as follows:
1. Generate B realizations of the gaussian $T \sim \mathcal{N}(0, \sigma^2)$, denoted by $T_1,..., T_B$
2.  Calculate $\sqrt{n}\left( \hat{F}\left(\hat{F}^{-1}(p) + \dfrac{T_b}{\sqrt{n}}\right) - p \right), b = 1,..., B$ and denote them as $y_b$, then the least squares estimate of $f(F^{-1}(p))$ is $\hat{A} = (x'x)^{-1}x'Y$, where $x= (T_1,..., T_B)^T$ and $Y = (y_1,..., y_B)^T$

### Variance selection 
#### Impact of the variance on the MSE
We advocate that the selection of the variance of the sampled Gaussians is crucial, as it can significantly affect the accuracy of the density estimation. 
In order to highlight this, we illustrate the estimation of the density at the median for an exponential distribution with rate $1.5$, and compute the MSE for a range of values of $\sigma$. We perform 100 repetitions of the code, assume censoring to be independent and to follow an exponential distribution with rate 0.48. We generate $B=10^5$ zero-mean Gaussians. 
The code to generate these simulations is available on `MSE/Densite_Point_Estim - Exponential.R`.

We notice that, especially for small sample sizes, the value of $\sigma$ plays an important role in the value of the MSE. Indeed, as sample size increases, the region where MSE is minimized becomes broader, which allows for greater flexibility when choosing the variance within a wider interval of values where the MSE is close to zero. For all sample sizes, we observe that the MSE has a pattern where it decreases until it reaches a plateau, where it remains low over a range of $\sigma$ values, before increasing again. Our goal in practical applications where survival and censoring distributions are unknown is to be able to automatically detect, from the observed times, such a plateau, and select a value of $\sigma$ that lays inside this interval, which grows broader as sample size increases.

![new_expo_n50_M100](https://github.com/user-attachments/assets/572cc906-e562-4f06-8847-dc4873499e58)

![new_expo_n200_M100](https://github.com/user-attachments/assets/7c72629f-3ff9-434e-a1b8-ca12dc407801)

![new_expo_n1000_M100](https://github.com/user-attachments/assets/ca0f7a68-b66a-423b-9898-464795b7f8b8)


#### Grid-search algorithm for variance selection
In practice the true value of the density at the quantile is unknown, so we cannot directly compute the MSE to identify the plateau described earlier. We propose a grid-search algorithm to automatically identify this plateau, which will allow us to select an appropriate value for $\sigma$. This procedure is implemented in `src/code_source` and `src/simulations_and_results`.

We illustrate this algorithm for an exponential distribution with rate 1.5 with censoring that follows an exponential with rate 0.12 (approximately 10% of censoring). We are interested in estimating the density of the exponential at the median, for which the true value is 0.75. We then perform the estimation of the density at the quantile for a grid of $\sigma$  that ranges from 0.1 to 10, in steps of length 0.05. We observe that there exists an interval of values of $\sigma$ such that the estimation is small before it and decreases after it as well. We aim to choose a value of $\sigma$ that lays inside this interval, because this is the interval that corresponds to the values close to the true value (represented by the dotted red line). Our code allows us to detect such interval and then it chooses a value of $\sigma$ that lays inside it. The chosen value is represented in the vertical blue line in the corresponding figures. In this illustrative example, we obtain an estimation of the density at the median equal to 0.745831, for a corresponding value of $\sigma$ equal to 3.75. 

![example of density estim at exponential](https://github.com/user-attachments/assets/51a4f61a-6b2e-4cbf-8f38-0a4ad79a8388)

![diff for density estim exponential](https://github.com/user-attachments/assets/2c19be49-ae6a-4f97-8454-845116ead13f)

The folder src has the corresponding codes for computing the estimation of the density for the Exponential and Cauchy distributions. The estimations using our proposed procedure are compared to the ones obtained using kernel density estimation (KDE) with bandwidth parameter approximated using leave-one-out cross-validation. Our implementation of LOOCV in the presence of censoring is available on `src/code_source.R`.

## Test of equality of quantiles
In our paper "Designing clinical trials for the comparison of single and multiple quantiles with right-censored data", we propose new power formulas for the comparison of one quantile between two treatment groups, as well as for the comparison of a collection of quantiles.
We have found that the variance of the test statistic depends on the estimation of the probability density function of the distribution of failure times at the quantile being tested.
In our paper, we illustrate the application of the test of equality of quantiles using the resampling-based procedure described above. 
The codes from this paper are made available in this repository.

## About the repository
This repository is organised as follows:
- `/src`:
    - code_source: auxiliary codes for test in presence of data and for planning clinical trial using the kernel density estimation and our resampling procedure, from our paper "Designing clinical trials for the comparison of single and multiple quantiles with right-censored data"
    - planning a clinical trial: codes for the section of Planning a clinical trial from our paper "Designing clinical trials for the comparison of single and multiple quantiles with right-censored data"
    - simulations_and_results: codes for the comparison of our method with KDE, detailed in our paper "A note on a reampling procedure for density estimation at quantiles"
    - test on OAK data: applying the density estimation for the univariate and multivariate test of equality of quantiles detailed in our paper "Designing clinical trials for the comparison of single and multiple quantiles with right-censored data"
 
- `/MSE`: codes for the figures for MSE varying with variance, for exponential distribution

# References
Farah, Beatriz, Aurélien Latouche, and Olivier Bouaziz. "A note on a resampling procedure for estimating the density at a given quantile." arXiv preprint arXiv:2509.02207 (2025).
Beatriz Farah, Olivier Bouaziz, Aurélien Latouche. Designing clinical trials for the comparison of single and multiple quantiles with right-censored data. 2025. ⟨hal-05052387v2⟩
Lin, C., Zhang, L., & Zhou, Y. (2015). Conditional quantile residual lifetime models for right censored data. Lifetime data analysis, 21, 75-96.
