require(survival)
library(dplyr)
library(survminer)
library(latex2exp)
library("lrstat")
library(hesim)
library(devtools)
library(mgcv)
library(zoo)
library(flexsurv)
source("~/GitHub/dens-estimation-at-quantile/src/code_source.R")
#rm(list = ls())
##############################################################################################################
# Simulations for paper density estimation

#############################################################################################################
# Exponential survival and exponential censoring, compared at the median

# Parameters

nis = c(50, 200)
#n0 = 200

shape0 = 1.0  # for exponential
scale0 = 1/1.5

rate = 0.48

rate_cens = c(1, 0.48, 0.12)

kernel = 'gaussian'

sd_lst = seq(0.1, 10, by = 0.05)
h_vals = seq(0.1, 0.5, by = 0.02)
B_sd = 1e+03

region_size = 10 # for n0 = 50
#region_size = 20 # for n0 = 200

medi0 = scale0*(log(2))^(1/shape0)
f_atmedi0 = (shape0/scale0)*(medi0/scale0)^(shape0-1)*exp(-(medi0/scale0)^shape0)

p = 0.5 # when we compare two distributions with the same parameters, we need to set a p (otherwise it will be NaN)
#q = 0.25 # for the kernel density estimation

ALPHA = 0.05

########################################################################################
set.seed(12)
nb_it = 500
all_fs_lin = data.frame()
all_fs_kde = data.frame()
all_sds_lin = data.frame()
all_hs_kde = data.frame()
it = 0

n0 = 200
rate = 0.12

while(it<nb_it){
  print(it)
  truetimes0_unsort = rweibull(n0,shape=shape0,scale=scale0) # True times: censored and uncensored (=time to event) 
  censtimes0_unsort = rexp(n0,rate=rate)#rep(Inf,n0)# # Censoring times
  #censtimes0_unsort = rep(Inf,n0) # without censoring
  
  status0_unsort = as.numeric(truetimes0_unsort <= censtimes0_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes0_unsort = pmin(truetimes0_unsort, censtimes0_unsort) # Observed times: min between true times and censored times
  
  DATA1 <- data.frame("x" = alltimes0_unsort, "event" = status0_unsort) #truetimes0_unsort
  Time1 <- alltimes0_unsort
  
  #lin_res = find_flat_regions(n0, DATA1, sd_lst, p, B_sd)
  
  out1 <- arrange(DATA1, x)
  
  alltimes1 = sort(alltimes0_unsort)
  status1 = out1$event[order(out1$x)]
  
  # Compute sigma values
  f_sigmas <- lapply(sd_lst, function(sd) {
    #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
    new_lin_estimation(DATA1, p = p, sd = sd, return_surv = FALSE, B = B_sd)
  })
  
  f_sigmas = unlist(f_sigmas)
  
  lin_res = find_local_maxima_with_region(f_sigmas, region_size)$ind_sigma[1] # maybe not do that, or verify
  if(length(lin_res) == 0){ # if there is no local maxima in f_sigmas, we look for a local minima in abs(diff(f_sigmas))
    #lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)
    lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  }
  
  lin_res1 = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  #lin_res = find_local_minima_with_region( abs(diff(f_sigmas, 2)) , region_size)$ind_sigma[1]
  
  index_lin = min(c(lin_res1, lin_res)) + 20 # sum 20 if the steps are by 0.05, sum 10 if they are of 0.1
  #maybe its actually the maximum instead of the minimum. see if we are being too conservative
  
  if(is.na(f_sigmas[index_lin])){
    break
  }
  
  xxi = find_quantile(DATA1)
  
  best_h = mise_kde_old(alltimes0_unsort, status0_unsort, h_vals, kernel = kernel)$optimal_h
  
  f_kde = our_kde(xxi, alltimes0_unsort, status0_unsort, best_h, kernel)
  
  all_fs_kde <-rbind(all_fs_kde, f_kde)
  all_fs_lin <-rbind(all_fs_lin, f_sigmas[index_lin])
  #all_fs_lin <-rbind(all_fs_lin, lin_res$f_sigmas[1])
  all_sds_lin <-rbind(all_sds_lin, sd_lst[index_lin])
  all_hs_kde <- rbind(all_hs_kde, best_h)
  
  it = it+1
}

names(all_fs_lin) = "f_estim_lin"
names(all_fs_kde) = "f_estim_kde"
names(all_sds_lin) = "sigma_lin"
names(all_hs_kde) = "h_kde"

bias_lin = mean(all_fs_lin$f_estim_lin) - f_atmedi0
var_lin = var(all_fs_lin$f_estim_lin)
MSE_lin <- bias_lin^2+var_lin

bias_kde = mean(all_fs_kde$f_estim_kde) - f_atmedi0
var_kde = var(all_fs_kde$f_estim_kde)
MSE_kde <- bias_kde^2+var_kde

results_df <- data.frame(
  Estimator = c("Lin", "KDE"),
  Bias = c(bias_lin, bias_kde),
  Variance = c(var_lin, var_kde),
  MSE = c(MSE_lin, MSE_kde)
)

results_df


#######################################################################################
# new test with the maximum instead of the minimum
set.seed(12)
nb_it = 500
all_fs_lin = data.frame()
all_fs_kde = data.frame()
all_sds_lin = data.frame()
all_hs_kde = data.frame()
it = 0

n0 = 200
rate = 0.12

while(it<nb_it){
  print(it)
  truetimes0_unsort = rweibull(n0,shape=shape0,scale=scale0) # True times: censored and uncensored (=time to event) 
  censtimes0_unsort = rexp(n0,rate=rate)#rep(Inf,n0)# # Censoring times
  #censtimes0_unsort = rep(Inf,n0) # without censoring
  
  status0_unsort = as.numeric(truetimes0_unsort <= censtimes0_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes0_unsort = pmin(truetimes0_unsort, censtimes0_unsort) # Observed times: min between true times and censored times
  
  DATA1 <- data.frame("x" = alltimes0_unsort, "event" = status0_unsort) #truetimes0_unsort
  Time1 <- alltimes0_unsort
  
  #lin_res = find_flat_regions(n0, DATA1, sd_lst, p, B_sd)
  
  out1 <- arrange(DATA1, x)
  
  alltimes1 = sort(alltimes0_unsort)
  status1 = out1$event[order(out1$x)]
  
  # Compute sigma values
  f_sigmas <- lapply(sd_lst, function(sd) {
    #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
    new_lin_estimation(DATA1, p = p, sd = sd, return_surv = FALSE, B = B_sd)
  })
  
  f_sigmas = unlist(f_sigmas)
  
  lin_res = find_local_maxima_with_region(f_sigmas, region_size)$ind_sigma[1] # maybe not do that, or verify
  if(length(lin_res) == 0){ # if there is no local maxima in f_sigmas, we look for a local minima in abs(diff(f_sigmas))
    #lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)
    lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  }
  
  lin_res1 = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  #lin_res = find_local_minima_with_region( abs(diff(f_sigmas, 2)) , region_size)$ind_sigma[1]
  
  index_lin = max(c(lin_res1, lin_res)) + 20 # sum 20 if the steps are by 0.05, sum 10 if they are of 0.1
  #maybe its actually the maximum instead of the minimum. see if we are being too conservative
  
  if(is.na(f_sigmas[index_lin])){
    break
  }
  
  xxi = find_quantile(DATA1)
  
  best_h = mise_kde_old(alltimes0_unsort, status0_unsort, h_vals, kernel = kernel)$optimal_h
  
  f_kde = our_kde(xxi, alltimes0_unsort, status0_unsort, best_h, kernel)
  
  all_fs_kde <-rbind(all_fs_kde, f_kde)
  all_fs_lin <-rbind(all_fs_lin, f_sigmas[index_lin])
  #all_fs_lin <-rbind(all_fs_lin, lin_res$f_sigmas[1])
  all_sds_lin <-rbind(all_sds_lin, sd_lst[index_lin])
  all_hs_kde <- rbind(all_hs_kde, best_h)
  
  it = it+1
}

names(all_fs_lin) = "f_estim_lin"
names(all_fs_kde) = "f_estim_kde"
names(all_sds_lin) = "sigma_lin"
names(all_hs_kde) = "h_kde"

bias_lin = mean(all_fs_lin$f_estim_lin) - f_atmedi0
var_lin = var(all_fs_lin$f_estim_lin)
MSE_lin <- bias_lin^2+var_lin

bias_kde = mean(all_fs_kde$f_estim_kde) - f_atmedi0
var_kde = var(all_fs_kde$f_estim_kde)
MSE_kde <- bias_kde^2+var_kde

results_df <- data.frame(
  Estimator = c("Lin", "KDE"),
  Bias = c(bias_lin, bias_kde),
  Variance = c(var_lin, var_kde),
  MSE = c(MSE_lin, MSE_kde)
)

results_df

plot(sd_lst, f_sigmas, xlab = 'Sigma', ylab = 'Estimation of the density at the median', main = "n = 200")
abline(v = sd_lst[index_lin], lty = 1, col = 'blue')
abline(h = f_atmedi0, lty = 2, col = 'red')
legend("topright", legend = c("Selected sigma by grid-search", "Real value of density at quantile"), 
       col = c("blue", "red"), lty = c(1, 2))

plot(sd_lst, abs(f_sigmas - f_atmedi0), xlab = 'Sigma', ylab = 'Absolute difference between real value and estimation', main = "n = 200")
abline(v = sd_lst[index_lin],lty = 1, col = 'blue')
legend("topright", legend = c("Selected sigma by grid-search"), 
       col = c("blue"), lty = c(1))
########################################################################################
# Cauchy survival with exponential censoring
# Heavy tailed distribution

# Compared at the median

# Parameters
nis = c(50, 200)


location_param <- 0  #x0 # For Cauchy
scale_param <- 1 #gamma

min_unif = 0 # For Uniform
max_unif = 3
rates = c(0.1, 0.7, 2.3) # 10% of censoring, 25% of censoring, 40% of censoring

n0 = 200
rate = 2.3

h = 0.2
kernel = 'gaussian'
sd_lst = seq(0.1, 10, by = 0.05)
h_vals = seq(0.1, 0.5, by = 0.02)
B_sd = 1e+03
tolerance = 0.01
#region_size = 10 # for n0 = 50
region_size = 20 # for n0 = 200

medi0 = location_param
f_atmedi0 = 1/(pi*scale_param*(1+( (medi0 - location_param)/scale_param)^2 )) 

p = 0.5 # when we compare two distributions with the same parameters, we need to set a p (otherwise it will be NaN)
#q = 0.25 # for the kernel density estimation

ALPHA = 0.05

########################################################################################
set.seed(12)
nb_it = 500
all_fs_lin = data.frame()
all_fs_kde = data.frame()
all_sds_lin = data.frame()
all_hs_kde = data.frame()
it = 0
while(it<nb_it){
  print(it)
  truetimes0_unsort = rcauchy(n0, location = location_param, scale = scale_param) # True times: censored and uncensored (=time to event) 
  censtimes0_unsort = rexp(n0,rate=rate)
  #censtimes0_unsort = runif(n0, min = min_unif, max = max_unif) #rep(Inf,n0)# # Censoring times
  #censtimes0_unsort = rep(Inf,n0) # without censoring
  
  status0_unsort = as.numeric(truetimes0_unsort <= censtimes0_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes0_unsort = pmin(truetimes0_unsort, censtimes0_unsort) # Observed times: min between true times and censored times
  
  DATA1 <- data.frame("x" = alltimes0_unsort, "event" = status0_unsort) #truetimes0_unsort
  Time1 <- alltimes0_unsort
  
  #lin_res = find_flat_regions(n0, DATA1, sd_lst, p, B_sd)
  out1 <- arrange(DATA1, x)
  
  alltimes1 = sort(alltimes0_unsort)
  status1 = out1$event[order(out1$x)]
  
  # Compute sigma values
  f_sigmas <- lapply(sd_lst, function(sd) {
    #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
    new_lin_estimation(DATA1, p = p, sd = sd, return_surv = FALSE, B = B_sd)
  })
  
  f_sigmas = unlist(f_sigmas)
  
  lin_res = find_local_maxima_with_region(f_sigmas, region_size)$ind_sigma[1] # maybe not do that, or verify
  if(length(lin_res) == 0){ # if there is no local maxima in f_sigmas, we look for a local minima in abs(diff(f_sigmas))
    #lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)
    lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  }
  
  lin_res1 = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
  #lin_res = find_local_minima_with_region( abs(diff(f_sigmas, 2)) , region_size)$ind_sigma[1]
  
  index_lin = min(c(lin_res1, lin_res)) + 20 # sum 20 if the steps are by 0.05, sum 10 if they are of 0.1
  #index_lin = lin_res$ind_sigma[1] + 20
  
  if(is.na(f_sigmas[index_lin])){
    break
  }
  
  xxi = find_quantile(DATA1)
  
  best_h = mise_kde_old(alltimes0_unsort, status0_unsort, h_vals, kernel = kernel)$optimal_h
  
  f_kde = our_kde(xxi, alltimes0_unsort, status0_unsort, best_h, kernel)
  
  all_fs_kde <-rbind(all_fs_kde, f_kde)
  all_fs_lin <-rbind(all_fs_lin, f_sigmas[index_lin])
  #all_fs_lin <-rbind(all_fs_lin, lin_res$f_sigmas[1])
  all_sds_lin <-rbind(all_sds_lin, sd_lst[index_lin])
  all_hs_kde <- rbind(all_hs_kde, best_h)
  
  
  it = it+1
}

names(all_fs_lin) = "f_estim_lin"
names(all_fs_kde) = "f_estim_kde"
names(all_sds_lin) = "sigma_lin"

bias_lin = mean(all_fs_lin$f_estim_lin) - f_atmedi0
var_lin = var(all_fs_lin$f_estim_lin)
MSE_lin <- bias_lin^2+var_lin

bias_kde = mean(all_fs_kde$f_estim_kde) - f_atmedi0
var_kde = var(all_fs_kde$f_estim_kde)
MSE_kde <- bias_kde^2+var_kde

results_df <- data.frame(
  Estimator = c("Lin", "KDE"),
  Bias = c(bias_lin, bias_kde),
  Variance = c(var_lin, var_kde),
  MSE = c(MSE_lin, MSE_kde)
)

results_df

#####################
# Visualization of true Cauchy survival: to show that we have heavy tail
x <- seq(0, 10, length.out = 1000)

# Compute the survival function (1 - CDF)
survival_function <- 1 - pcauchy(x, location_param, scale_param)

# Plot the survival function
plot(x, survival_function, type = "l", col = "blue", lwd = 2,
     main = "Survival Function of Cauchy Distribution",
     xlab = "x", ylab = "Survival Probability S(x)", ylim = c(-0.2, 0.6))
normal_survival <- 1 - pnorm(x, mean = 0, sd = 1)
lines(x, normal_survival, col = "red", lwd = 2, lty = 2)

legend("topright", legend = c("Cauchy", "Normal"), col = c("blue", "red"),
       lty = c(1, 2), lwd = 2)
