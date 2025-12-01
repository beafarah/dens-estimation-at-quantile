# Resampling method for density estimation at quantile and Clinical trial design for the comparison of single and multiple quantiles with right-censored data
This repository contains functions and code developped for the following papers:

* *Designing clinical trials for the comparison of single and multiple quantiles with right-censored data*
* *A note on a resampling procedure for density estimation at quantiles*

## Overview:
### 📁 MSE

This folder includes:
- **`Densite_Point_Estim - Exponential.R`**: A script containing the simulation code for evaluating the Mean Squared Error (MSE) of density point estimation under an exponential distribution for various values of ( \sigma ).
- **`n50/`** ,  **`n200/`** ,  **`n1000/`** : Folders that store the outputs of the MSE simulations.

### 📁 scr

This folder contains all core scripts for resampling, testing, and simulations.

- **`code_source.R`**: Script with all functions necessary for the resampling procedure for density estimation, as well as for the univariate and multivariate tests of equality of quantiles, with the explicit power from the analytical formula for designing clinical trials and the estimated power for the test in the presence of data.
  Some of the main functions:
```
- SurvQtestExp_new & SurvQtestPieceExp_new: analytical power for the univariate test for a given difference in quantiles Δ for exponential and piecewise exponential treatments, respectively
- SurvQtestExp_empirical & SurvQtestPieceExp_empirical: empirical rejection rate for exponential and piecewise exponential treatments, respectively
- new_lin_estimation: estimation of the density at a fixed quantile for a given value of σ
```

- **`planning a clinical trial.R`**: Includes all simulation code used in the *Planning a clinical trial* section of the paper *Designing clinical trials for the comparison of single and multiple quantiles with right-censored data*
-  **`test on OAK data.R`**: Provides the results of applying the univariate and multivariate tests to data from the OAK clinical trial, which correspond to section *Application of the test on data from the OAK study* of the paper *Designing clinical trials for the comparison of single and multiple quantiles with right-censored data*
-  **`simulations_and_results.R`**: Contains all simulations related to the paper *A note on a resampling procedure for density estimation at quantiles*, which correspond to the resampling method to estimate the density at a given quantile.

<!--
## Resampling procedure
Let $\hat{F}$ denote the Kaplan-Meier estimator of the distribution $F$, and $\hat{F}^{-1}(p)$ its inverse at probability $p$. The procedure to estimate $f(\hat{F}^{-1}(p))$ is as follows:
1. Generate B realizations of the gaussian $T \sim \mathcal{N}(0, \sigma^2)$, denoted by $T_1,..., T_B$
2.  Calculate $\sqrt{n}\left( \hat{F}\left(\hat{F}^{-1}(p) + \dfrac{T_b}{\sqrt{n}}\right) - p \right), b = 1,..., B$ and denote them as $y_b$, then the least squares estimate of $f(F^{-1}(p))$ is $\hat{A} = (x'x)^{-1}x'Y$, where $x= (T_1,..., T_B)^T$ and $Y = (y_1,..., y_B)^T$

### Variance selection 
#### Impact of the variance on the MSE
We advocate that the selection of the variance of the sampled Gaussians is crucial, as it can significantly affect the accuracy of the density estimation. 
In order to highlight this, we illustrate the estimation of the density at the median for an exponential distribution with rate $1.5$, and compute the MSE for a range of values of $\sigma$. We perform 100 repetitions of the code, assume censoring to be independent and to follow an exponential distribution with rate 0.48. We generate $B=10^5$ zero-mean Gaussians. 
The code to generate these simulations is available on `MSE/Densite_Point_Estim - Exponential.R`.

We notice that, especially for small sample sizes, the value of $\sigma$ plays an important role in the value of the MSE. Indeed, as sample size increases, the region where MSE is minimized becomes broader, which allows for greater flexibility when choosing the variance within a wider interval of values where the MSE is close to zero. For all sample sizes, we observe that the MSE has a pattern where it decreases until it reaches a plateau, where it remains low over a range of $\sigma$ values, before increasing again. Our goal in practical applications where survival and censoring distributions are unknown is to be able to automatically detect, from the observed times, such a plateau, and select a value of $\sigma$ that lays inside this interval, which grows broader as sample size increases.

<p float="left">
  <img src="https://github.com/user-attachments/assets/572cc906-e562-4f06-8847-dc4873499e58" alt="n50" width="300" style="margin-right:10px;">
  <img src="https://github.com/user-attachments/assets/7c72629f-3ff9-434e-a1b8-ca12dc407801" alt="n200" width="300" style="margin-right:10px;">
  <img src="https://github.com/user-attachments/assets/ca0f7a68-b66a-423b-9898-464795b7f8b8" alt="n1000" width="300">
</p>


#### Grid-search algorithm for variance selection
In practice the true value of the density at the quantile is unknown, so we cannot directly compute the MSE to identify the plateau described earlier. We propose a grid-search algorithm to automatically identify this plateau, which will allow us to select an appropriate value for $\sigma$. This procedure is implemented in `src/code_source` and `src/simulations_and_results`.

We illustrate this algorithm for an exponential distribution with rate 1.5 with censoring that follows an exponential with rate 0.12 (approximately 10% of censoring). We are interested in estimating the density of the exponential at the median, for which the true value is 0.75. We then perform the estimation of the density at the quantile for a grid of $\sigma$  that ranges from 0.1 to 10, in steps of length 0.05. We observe that there exists an interval of values of $\sigma$ such that the estimation is small before it and decreases after it as well. We aim to choose a value of $\sigma$ that lays inside this interval, because this is the interval that corresponds to the values close to the true value (represented by the dotted red line). Our code allows us to detect such interval and then it chooses a value of $\sigma$ that lays inside it. The chosen value is represented in the vertical blue line in the corresponding figures. In this illustrative example, we obtain an estimation of the density at the median equal to 0.745831, for a corresponding value of $\sigma$ equal to 3.75. 

<p float="left">
  <img src="https://github.com/user-attachments/assets/51a4f61a-6b2e-4cbf-8f38-0a4ad79a8388" alt="density estimation example" width="400" style="margin-right:10px;">
  <img src="https://github.com/user-attachments/assets/2c19be49-ae6a-4f97-8454-845116ead13f" alt="density estimation difference" width="400">
</p>

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


-->

# References
Farah, Beatriz, Aurélien Latouche, and Olivier Bouaziz. "A note on a resampling procedure for estimating the density at a given quantile." arXiv preprint arXiv:2509.02207 (2025).

Beatriz Farah, Olivier Bouaziz, Aurélien Latouche. Designing clinical trials for the comparison of single and multiple quantiles with right-censored data. 2025. ⟨hal-05052387v2⟩

Lin, C., Zhang, L., & Zhou, Y. (2015). Conditional quantile residual lifetime models for right censored data. Lifetime data analysis, 21, 75-96.
