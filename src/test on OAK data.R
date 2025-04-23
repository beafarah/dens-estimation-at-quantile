library(survival)
library(survminer)
setwd("C:/Users/beafa/OneDrive/Documents/GitHub/dens-estimation-at-quantile/src")
load("data_ICI_Rittmeyer.RData")
#rm(list = ls())
#View(data_ICI_Rittmeyer)
# event = 1 for death and 0 for censoring
#source("~/GitHub/these/kosorok/kosorok_test_source.R")
source("code_source.R")

########################################################
## Codes for test on OAK data

########################################################
# Univariate test

fit_KM = survfit(Surv(time, event)~tmt.arm.number, data = data_ICI_Rittmeyer)

ggsurvplot(fit_KM, risk.table = TRUE, xlim = c(0.4, 27), break.time.by = 3, legend.labs = c("Chemo", "Immuno"))

########################################################
data_chemo = data_ICI_Rittmeyer[data_ICI_Rittmeyer$tmt.arm.number == 0,]
data_immuno = data_ICI_Rittmeyer[data_ICI_Rittmeyer$tmt.arm.number == 1,]
# 89 censored for chemo,
# 206 censored for immuno

##########################################################
data_group0 = select(data_chemo, time, event)
names(data_group0) = c("x", "event")
row.names(data_group0) = NULL
data_group1 = select(data_immuno, time, event)
names(data_group1) = c("x", "event")
row.names(data_group1) = NULL

Time0 = data_group0$x
Time1 = data_group1$x

alltimes0_sort = sort(data_group0$x)
status_sort0 = data_group0$event[order(data_group0$x)]

alltimes = alltimes0_sort
status = status_sort0

uncensored_times <- alltimes[status == 1]

###########
# leave-one-out cross validation function
h_minimizer = function(h){
  "function that will be minimized for the cross validation. returns one value for one bandwidth h."
  num <- length(alltimes)
  G_hat <- survfit(Surv(alltimes, 1 - status) ~ 1) # KM estimator for censoring survival
  scalier_fun <- stepfun(G_hat$time, c(1, G_hat$surv))
  
  aux_fun = function(t){
    soma = 0
    for(uncens_t in uncensored_times){
      soma = soma + K_gaussian((uncens_t - t) / h) / scalier_fun(uncens_t)
    }
    return( (1/(num*h)* soma)^2 )
  } 
  bounds = c(min(alltimes), max(alltimes))
  first_term = integrate(aux_fun, lower = min(alltimes), upper = max(alltimes))
  
  for(i in 1:(length(uncensored_times) - 1)){
    som = 0 
    for(j in 1:(length(uncensored_times) - 1)){
      if(i != j){
        som = som + K_gaussian( (uncensored_times[i] - uncensored_times[j])/h ) /( scalier_fun(uncensored_times[i]) * scalier_fun(uncensored_times[j]))
      }
    }
    som = -2/(num*(num-1)*h)*som
  }
  
  return(first_term$value + som)
  
}

#####################################################################
# choice of bandwidth that minimizes the function for cross validation
#h_vals = seq(0.1, 0.5, by = 0.02)
h_vals = seq(0.1, 1, by = 0.1)
mises = c()
for(h in h_vals){
  mises = c(mises, h_minimizer(h))
}

index_h_opt = which.min(mises)
h_opt = h_vals[index_h_opt]

##########################################################################
# pvalue for kde and resampling
p = 0.5

kernel = 'gaussian'
sd_lst = seq(0.1, 10, by = 0.05)
h_vals = seq(0.1, 0.5, by = 0.02)
B = 1e+03
B_sd = 1e+03
tolerance = 0.01
region_size = 10

# grid search for sigma selection (resampling method) and cross validation for bandwidth selection (kde)
f_sigmas <- lapply(sd_lst, function(sd) {
  #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
  new_lin_estimation(data_group0, p = p, sd = sd, return_surv = FALSE, B = B_sd)
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
sd_auto = sd_lst[index_lin]

####################
sim_lin_auto = SurvQtestLin(data_group0, Time0, data_group1, Time1,p=p, sd = sd_auto, B = B)
sim_kde_auto = SurvQtest_kde(data_group0, Time0, data_group1, Time1,p=p, h = h_opt, ker = kernel)


###################
# compare with Tang and Jeong at the median (application of Brookmeyer)
res_TJ = pval_Tang_Jeong(data_group0, Time0, data_group1, Time1, p, return_stat = FALSE)


#################################################
# Multivariate test
########################################################
# Estimation of the difference of quantiles
p1 = 0.7
p2 = 0.5
alltimes1_unsort = data_group0$x
status1_unsort = data_group0$event
alltimes2_unsort = data_group1$x
status2_unsort = data_group1$event

km_fit_F <- survfit(Surv(alltimes1_unsort, status1_unsort) ~ 1)
quantiles_F <- quantile(km_fit_F, probs = c(p1, p2))$quantile

km_fit_G <- survfit(Surv(alltimes2_unsort, status2_unsort) ~ 1)
quantiles_G <- quantile(km_fit_G, probs = c(p1, p2))$quantile

deltas_hat = quantiles_F - quantiles_G
print(deltas_hat)


#######################################################
# choice of sd for lin
p1 = 0.7
p2 = 0.5

kernel = 'gaussian'
sd_lst = seq(0.1, 10, by = 0.05)
h_vals = seq(0.1, 0.5, by = 0.02)
B = 1e+03
B_sd = 1e+03
tolerance = 0.01
region_size = 10

##
# sd for p1
f_sigmas <- lapply(sd_lst, function(sd) {
  #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
  new_lin_estimation(data_group0, p1, sd = sd, return_surv = FALSE, B = B_sd)
})

f_sigmas = unlist(f_sigmas)

lin_res = find_local_minima_with_region(f_sigmas, region_size)$ind_sigma[1] # maybe not do that, or verify
if(length(lin_res) == 0){ # if there is no local maxima in f_sigmas, we look for a local minima in abs(diff(f_sigmas))
  #lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)
  lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
}

lin_res1 = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
#lin_res = find_local_minima_with_region( abs(diff(f_sigmas, 2)) , region_size)$ind_sigma[1]

index_lin = min(c(lin_res1, lin_res)) + 20 # sum 20 if the steps are by 0.05, sum 10 if they are of 0.1
sd_auto1 = sd_lst[index_lin]

##
# sd for p2
f_sigmas <- lapply(sd_lst, function(sd) {
  #lin_estimation(out1, alltimes1, p, sd = sd, return_surv = FALSE, B = B_sd)
  new_lin_estimation(data_group0, p = p2, sd = sd, return_surv = FALSE, B = B_sd)
})

f_sigmas = unlist(f_sigmas)

lin_res = find_local_minima_with_region(f_sigmas, region_size)$ind_sigma[1] # maybe not do that, or verify
if(length(lin_res) == 0){ # if there is no local maxima in f_sigmas, we look for a local minima in abs(diff(f_sigmas))
  #lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)
  lin_res = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
}

lin_res1 = find_local_minima_with_region( abs(diff(f_sigmas)) , region_size)$ind_sigma[1]
#lin_res = find_local_minima_with_region( abs(diff(f_sigmas, 2)) , region_size)$ind_sigma[1]

index_lin = min(c(lin_res1, lin_res)) + 20 # sum 20 if the steps are by 0.05, sum 10 if they are of 0.1
sd_auto2 = sd_lst[index_lin]

SurvQtestLin_Multi(data_group0, Time0, data_group1, Time1, p1, p2, sd = sd_auto1)

## kde
###########
# choice of h for kde
alltimes0_sort = sort(data_group0$x)
status_sort0 = data_group0$event[order(data_group0$x)]
#best_h = lscv_optimal_bandwidth(alltimes = alltimes0_sort, status = status_sort0, bandwidth_range = h_vals, kernel = kernel)
#best_h = mise_kde_old(alltimes0_sort, status_sort0, h_vals, kernel = kernel)$optimal_h

h_vals = seq(0.1, 1, by = 0.1)
mises = c()
for(h in h_vals){
  mises = c(mises, h_minimizer(h))
}

index_h_opt = which.min(mises)
h_opt = h_vals[index_h_opt]
######################
#p1 = 0.0001
p1 = 0.7
#p2 = 0.0005 #0.5, 0.75
p2 = 0.5


SurvQtestLin_Multi_kde(data_group0, Time0, data_group1, Time1, p1, p2, h = h_opt)


surv_plot = ggsurvplot(fit_KM, risk.table = TRUE, xlim = c(0.4, 27), break.time.by = 3, legend.labs = c("Chemo", "Immuno"))


surv_plot$plot <- surv_plot$plot + labs(color = "Treatment", linetype = "Treatment", shape = "Treatment")+
  geom_hline(yintercept = 1-0.1, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 1-0.05, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 1-0.5, linetype = "dashed", color = "grey")+
  geom_hline(yintercept = 1-0.7, linetype = "dashed", color = "grey")

surv_plot$table <- surv_plot$table +
  ylab("Treatment") 

print(surv_plot)

