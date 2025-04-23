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
#rm(list = ls())
###############################
# auxiliary functions
g_piecewise = function(tcut, lamb_begin, lamb_end, x){
  if(x < tcut){
    res = lamb_begin*exp(-lamb_begin*x)
  }else{
    res = lamb_end*exp(-lamb_end*(x - tcut) - lamb_begin*tcut)
  }
  return(res)
}

S_piecewise = function(tcut, lamb_begin, lamb_end, x){
  if(x < tcut){
    res = exp(-lamb_begin*x)
  }else{
    res = exp(-lamb_begin*tcut)*exp(-lamb_end*(x-tcut))
  }
  return(res)
}


data_proportional_piece_Exp <- function(n1, n2, k1, k2, lamb_group1, tcut, p, delta, late = TRUE, cens_type = 'unif', cens_rate = 0.3){
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  data1 <- data.frame("x" = alltimes1_unsort , "event" = status1_unsort) #truetimes0_unsort
  Time1 <- alltimes1_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  quantile_group1 = -log(1-p)/lamb_group1
  quantile_group2 = -delta + quantile_group1
  
  if(late){
    lamb_begin = lamb_group1
    if(p>=1-exp(-lamb_begin*tcut) & tcut-quantile_group2!=0){
      lamb_group2 = (log(1-p)+lamb_begin*tcut)/(tcut-quantile_group2)
      if(lamb_group2<0){
        lamb_group2 = -log(1-p)/quantile_group2
      }
    }
    else{
      print('warning: values of p and tcut are not compatible') 
      lamb_group2 = -log(1-p)/quantile_group2
      if(lamb_group2<0){
        lamb_group2 = (log(1-p)+lamb_begin*tcut)/(tcut-quantile_group2)
      }
    }
    lamb_end = lamb_group2
  }
  else{
    lamb_end = lamb_group1
    lamb_begin_candidate1 = -log(1-p)/quantile_group2
    lamb_begin_candidate2 = (lamb_end*(tcut - quantile_group2) - log(1-p))/tcut
    if(p<1-exp(-lamb_begin_candidate1*tcut)){
      lamb_group2 = lamb_begin_candidate1
    }
    else{
      lamb_group2 = lamb_begin_candidate2
    }
    lamb_begin = lamb_group2
  }
  
  
  
  #print(c('lamb_begin' = lamb_begin, 'lamb_end' = lamb_end))
  
  rate = c(lamb_begin, lamb_end)
  
  truetimes2_unsort = rsurv(n2, cuts = c(tcut), alpha = rate)
  
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n2, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n2, rate = cens_rate)
  }
  
  #censtimes2_unsort <- runif(n2, min = k2, max = k1+k2)
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  data2 <- data.frame("x" = alltimes2_unsort, "event" = status2_unsort)
  Time2 <- alltimes2_unsort
  
  data2_PLOT = data.frame(data2)
  data2_PLOT$groups = rep(2, n2)
  
  DATA_PLOT = rbind(data1_PLOT, data2_PLOT)
  
  
  return(list("data1" = data1, "Time1" = Time1, "data2"= data2, "Time2" = Time2, "DATA_PLOT" = DATA_PLOT,
              "quantile_group1" = quantile_group1, "quantile_group2" = quantile_group2, "lamb_group2" = lamb_group2))
  
}

###############################
# analytical power

# exponential treatment
analyticalpwr_exp = function(n1, n2, p, lamb_group1, delta, cens_rate, alpha = 0.05){
  n = n1 + n2
  mu = n1/n
  quantile_group1 = -log(1-p)/lamb_group1
  quantile_group2 = -delta + quantile_group1
  lamb_group2 = -log(1-p)/quantile_group2
  f = lamb_group1*exp(-lamb_group1 * quantile_group1)
  g = lamb_group2*exp(-lamb_group2 * quantile_group2)
  
  phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
  gamma = (1-p)^2*(lamb_group2/(lamb_group2+cens_rate))*(exp((lamb_group2+cens_rate)*quantile_group2) -1)
  
  var_H0 = phi/(mu*(f^2))  + gamma/((1-mu)* g^2)
  
  z <- qnorm(1 - alpha / 2)
  
  pwr = 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  return(pwr)
}

# piecewise treatment
analyticalpwr_piece = function(n1, n2, p, lamb_group1, tcut, delta, cens_rate, alpha = 0.05){
  n = n1 + n2
  mu = n1/n
  quantile_group1 = -log(1-p)/lamb_group1
  quantile_group2 = -delta + quantile_group1
  if(p > 1-exp(-lamb_group1*tcut)){
    lamb_end = (log(1-p) + lamb_group1*tcut)/(delta - quantile_group1 + tcut)
  }
  
  else{
    print('verify values for p and tcut')
  }
  
  f = lamb_group1*exp(-lamb_group1 * quantile_group1)
  #g = lamb_group2*exp(-lamb_group2 * quantile_group2)
  
  g = g_piecewise(tcut, lamb_group1, lamb_end, quantile_group2)
  
  phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
  gamma = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*tcut) -1) +( ((1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_group1-lamb_end)*tcut)) )*((exp((lamb_end+cens_rate)*(quantile_group2))) - exp((lamb_end+cens_rate)*(tcut))) )
  
  var_H0 = phi/(mu*(f^2))  + gamma/((1-mu)* g^2)
  
  z <- qnorm(1 - alpha / 2)
  
  pwr = 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  return(pwr)
}

##################################################
# 1) Planning a clinical trial

###############################################
# exponential control, piecewise exponential treatment
SurvQtestPieceExp_empirical <- function(n1, n2, k1, k2, lamb_group1, lamb_end, tcut, p, cens_type = 'unif', cens_rate = 0.3, alpha = 0.05, late = TRUE) {
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  #####
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  data1 <- data.frame("x" = alltimes1_unsort , "event" = status1_unsort) #truetimes0_unsort
  #Time1 <- alltimes1_unsort
  data1 = data1[order(data1$x), ]
  Time1 = data1$x
  ##
  
  #status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  #alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  #Tsort1 = sort(alltimes1_unsort, index.return = TRUE)
  #alltimes1 = Tsort1$x
  #status1 = alltimes1[Tsort1$ix]
  #data1 <- data.frame("x" = alltimes1, "event" = status1) #truetimes0_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  #quantile_group1 = find_quantile(data1)
  km_fit1 <- survfit(Surv(alltimes1_unsort, status1_unsort) ~ 1)
  quantile_group1 = unname(summary(km_fit1)$table["median"])
  
  lamb_begin = lamb_group1
  
  rate = c(lamb_begin, lamb_end)
  truetimes2_unsort = rsurv(n2, cuts = c(tcut), alpha = rate)
  
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n2, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n2, rate = cens_rate)
  }
  
  #status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  #alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  #Tsort2 = sort(alltimes2_unsort, index.return = TRUE)
  #alltimes2 = Tsort2$x
  #status2 = alltimes2[Tsort2$ix]
  
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  data2 <- data.frame("x" = alltimes2_unsort , "event" = status2_unsort) #truetimes0_unsort
  #Time2 <- alltimes2_unsort
  data2 = data2[order(data2$x), ]
  Time2 = data2$x
  ##
  
  #data2 <- data.frame("x" = alltimes2, "event" = status2) #truetimes0_unsort
  
  #Time2 <- alltimes2_unsort
  
  km_fit2 <- survfit(Surv(alltimes2_unsort, status2_unsort) ~ 1)
  quantile_group2 = unname(summary(km_fit2)$table["median"])
  
  #quantile_group2 = find_quantile(data2)
  delta = quantile_group1 - quantile_group2
  
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  
  #fun_integrand = function(rate, x){
  #return(rate*x*exp(rate*x))
  #}
  
  # estimating phi and gamma numerically
  #phi = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group1, rate = lamb_group1)$val # they are the same as computing with the 'true_phi' function
  #gamma = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group2, rate = lamb_group2)$val
  
  # computing explicit phi and gamma
  if(cens_type == 'unif'){
    if(quantile_group1 <k2){
      phi = (1-p)^2*(exp(lamb_group1*quantile_group1) - 1)
      cat('ERROR: Add computation for gamma')
      stop("ERROR")
    }
  }
  
  else if(cens_type == 'exp'){
    phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
    #gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) + (1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin+cens_rate)*tcut))*(exp((lamb_end+cens_rate)*(quantile_group2 - tcut)))
    gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) +( ((1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin-lamb_end)*tcut)) )*((exp((lamb_end+cens_rate)*(quantile_group2))) - exp((lamb_end+cens_rate)*(tcut))) )
    
  }
  
  
  f = lamb_group1*exp(-lamb_group1 * quantile_group1)
  g = g_piecewise(tcut, lamb_begin, lamb_end, quantile_group2)
  
  var_H0 = phi/(mu*(f^2))  + gamma/((1-mu)* g^2)
  
  stat = (delta)/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  #PVALUE <- 1 - pchisq(stat^2, df = 1) # because the stat is squared, we use chisq instead of gaussian
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = f, 'g' = g, 'gamma' = gamma, 'delta' = delta))
  
  result <- data.frame(
    NUM_F = n1,
    NUM_G = n2,
    XI_F = quantile_group1,
    XI_G = quantile_group2,
    delta = delta,
    p = p,
    LAMB_GROUP2 = lamb_end,
    STAT = stat,
    g = g,
    GAMMA = gamma,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    PVALUE = PVALUE, 
    POWER = POWER
  )
  
  return(result)
  
}
##################################################################
# sample size calculation

# exponential with censoring:
samplesize_Exp_new = function(lamb_group1, p, delta,cens_rate, n_initial1 = 10, n_initial2 = 10, n_step = 1, alpha = 0.05, target_power = 0.95, n_lim = 5000, return_list = FALSE){
  n_group1 = n_initial1
  n_group2 = n_initial2
  pwr_estim = analyticalpwr_exp(n_group1, n_group2, p, lamb_group1, delta, cens_rate)
  all_pwrs = c(pwr_estim)
  all_ns = c(n_group1)
  while(n_group1 < n_lim && n_group2 <n_lim && pwr_estim < target_power){
    n_group1 = n_group1 + n_step
    n_group2 = n_group2 + n_step
    pwr_estim = analyticalpwr_exp(n_group1, n_group2, p, lamb_group1, delta, cens_rate)
    all_pwrs = c(all_pwrs, pwr_estim)
    all_ns = c(all_ns, n_group1)
  }
  if(n_group1 == n_lim || n_group2 == n_lim){
    print('maximum sample size attained')
  }
  if(return_list){
    return(list(pwr_estim = pwr_estim, n_group1 = n_group1, all_pwrs = all_pwrs, all_ns = all_ns) )
  }
  return(c(pwr_estim = pwr_estim, n_group1 = n_group1))
}

# piecewise exponential with censoring
samplesize_PieceExp_new = function(lamb_group1, tcut, p, delta, cens_rate, n_initial1 = 10, n_initial2 = 10, n_step = 1, alpha = 0.05, target_power = 0.95, n_lim = 5000, return_list = FALSE){
  n_group1 = n_initial1
  n_group2 = n_initial2
  pwr_estim = analyticalpwr_piece(n_group1, n_group2, p, lamb_group1, tcut, delta, cens_rate)
  all_pwrs = c(pwr_estim)
  all_ns = c(n_group1)
  while(n_group1 < n_lim && n_group2 <n_lim && pwr_estim < target_power){
    n_group1 = n_group1 + n_step
    n_group2 = n_group2 + n_step
    pwr_estim = analyticalpwr_piece(n_group1, n_group2, p, lamb_group1, tcut, delta, cens_rate)
    all_pwrs = c(all_pwrs, pwr_estim)
    all_ns = c(all_ns, n_group1)
  }
  if(n_group1 == n_lim || n_group2 == n_lim){
    print('maximum sample size attained')
  }
  if(return_list){
    return(list(pwr_estim = pwr_estim, n_group1 = n_group1, all_pwrs = all_pwrs, all_ns = all_ns) )
  }
  return(c(pwr_estim = pwr_estim, n_group1 = n_group1))
}



#####################################################################
# Test in the presence of data

# lin resampling estimation
new_lin_estimation <- function(data1, p, sd = 1, return_surv = FALSE, B = 100) {
  # function that computes the estimation of density at the quantiles using the method from Lin (2015).
  # note that data1 should be sorted from the oldest to most recent time!!!!
  # INPUT:
  # data1 (Required) = Random sample from F. It consists of two columns. First column contains the observed times. Second column contains censoring value( 1 = event, 0 = censoring ).
  # Time1 (Required) = Variable containing time to event or last follow-up in any units in DATA1.
  # p = The probability of which the quantile would be tested. The null hypothesis is that the pth quantile of the F distribution is same as the pth quantile of the G distribution. .
  # sd = Standard deviation for the normal distribution (default is 1)
  # B = number of normal that are generated
  # OUTPUT:
  # if return_surv: list with density estimation at the quantile and the KM estimation for the survival function
  # otherwise, returns the density estimation at the quantile
  
  ##########
  xxi = find_quantile(data1)
  
  # Kaplan Meier estimator of the survival
  F_hat = survfit(Surv(x, event == 1) ~ 1, data = data1) # KM estimator of the survival 
  F_scalier_fun <- stepfun(F_hat$time,c(1,F_hat$surv))
  
  # distribution function
  F_dist = function(x){
    1- F_scalier_fun(x) 
  }
  
  ########################
  n0 = nrow(data1)
  Zi_vec1 = rnorm(n = B, mean = 0, sd = sd)
  Yi_vec1 = sqrt(n0)*( F_dist( xxi + Zi_vec1/sqrt(n0) ) - p )
  model = lm( Yi_vec1 ~ (Zi_vec1) -1 )
  
  f = unname(coef(model))
  
  if(return_surv == TRUE){
    return(c(F_scalier_fun, model) ) # return the survival function and the estimation of survival function 
  }
  return(f)
}


#########################################################
# grid search for choosing the variance of the generated gaussians
find_local_maxima_with_region <- function(values, region_size = 1) {
  # Check if the input vector has enough elements for the given region
  if (length(values) < (2 * region_size + 1)) {
    stop("The input vector is too short for the specified region size.")
  }
  
  # Initialize empty lists to store indices and values of local maxima
  maxima_indices <- c()
  maxima_values <- c()
  
  # Loop through the vector, considering only valid indices for the region
  for (i in (region_size + 1):(length(values) - region_size)) {
    # Extract the region around the current element
    region <- values[(i - region_size):(i + region_size)]
    
    # Check if the current element is the maximum in the region
    if (values[i] == max(region) && sum(values[i] == region) == 1) {
      maxima_indices <- c(maxima_indices, i)
      maxima_values <- c(maxima_values, values[i])
    }
  }
  return(data.frame(ind_sigma = maxima_indices, f_sigmas = maxima_values))
}

# region_size determines the size of the neighborhood around each point that we examine to check if it's a local maximum and if there's a plateau around it.

find_local_minima_with_region <- function(values, region_size = 1) {
  # Check if the input vector has enough elements for the specified region
  if (length(values) < (2 * region_size + 1)) {
    stop("The input vector is too short for the specified region size.")
  }
  
  # Initialize empty lists to store indices and values of local minima
  minima_indices <- c()
  minima_values <- c()
  
  # Loop through the vector, considering only valid indices for the region
  for (i in (region_size + 1):(length(values) - region_size)) {
    # Extract the region around the current element
    region <- values[(i - region_size):(i + region_size)]
    
    # Check if the current element is the minimum in the region
    if (values[i] == min(region) && sum(values[i] == region) == 1) {
      minima_indices <- c(minima_indices, i)
      minima_values <- c(minima_values, values[i])
    }
  }
  
  # Return a data frame with indices and values of local minima
  return(data.frame(ind_sigma = minima_indices, abs_diff_f = minima_values))
}

find_quantile = function(DATA1){
  n0 = nrow(DATA1)
  Time1 = DATA1$x
  out1 <- arrange(DATA1, Time1)
  x <- as.matrix(out1)
  xtmp1 <- matrix(0, n0, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n0) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n0)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n0) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n0)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n0) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  xtmp4 <- numeric(n0)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n0) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n0)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n0) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  # single numeric value representing the minimum index in the xtmp4 vector where the value is greater than or equal to the specified threshold p
  index1 <- min(which(xtmp4 >= p))
  
  # value from the first column of the matrix x at the row index1
  xxi <- x[index1, 1]
  return(unname(xxi))
}


################################################
# KDE
# Using kernel density estimation
K_triang <- function(x){
  # triangular ker 
  y=x
  y[-1<=x & x<=0] <- x[-1<=x & x<=0]+1
  y[x<=1 & x>0] <- 1-x[x<=1 & x>0]
  y[x<(-1) | x>1] <- rep(0,sum(x<(-1) | x>1))
  return(y)
}

K_gaussian <- function(x){
  # gaussian ker
  return((1/sqrt(2*pi))*exp(-(1/2)*x^2))
}

our_kde<-function(t, alltimes, status, h, ker = "triangular"){
  # alltimes = observed (censored) times sorted
  # status = censoring status sorted
  num = length(alltimes)
  G_hat = survfit(Surv(alltimes, 1- status)~1) # KM estimator of the survival function of censoring 
  scalier_fun <- stepfun(G_hat$time,c(1,G_hat$surv))
  uncensored_times = alltimes[status == 1]
  
  if(ker == "triangular"){
    kde = sum(K_triang((uncensored_times-t)/h)/(scalier_fun(uncensored_times)))/(num*h)
  }
  else if(ker == "gaussian"){
    kde = sum(K_gaussian((uncensored_times-t)/h)/(scalier_fun(uncensored_times)))/(num*h)
  }
  else{
    return("ker should be triangular or gaussian")
  }
  return(kde)
}

# cross validation
lscv_optimal_bandwidth <- function(alltimes, status, bandwidth_range, kernel = "triangular") {
  # leave one out cv
  n <- length(alltimes)
  
  # Calculate LSCV values for each bandwidth in the range
  lscv_values <- sapply(bandwidth_range, function(bandwidth) {
    lscv_sum <- 0
    
    # Sum up leave-one-out KDE estimates
    for (i in 1:n) {
      alltimes_loo <- alltimes[-i]
      status_loo <- status[-i]
      
      kde_estimate <- our_kde(alltimes[i], alltimes_loo, status_loo, h = bandwidth, ker = kernel)
      lscv_sum <- lscv_sum + kde_estimate
    }
    
    # Compute the full KDE for the integrated squared error part
    kde_full <- sapply(alltimes, function(t) our_kde(t, alltimes, status, h = bandwidth, ker = kernel))
    l2_norm <- mean(kde_full^2)
    
    # Calculate and return the LSCV criterion value for this bandwidth
    lscv_value <- l2_norm - (2 / n) * lscv_sum
    return(lscv_value)
  })
  
  # Find the optimal bandwidth with the minimum LSCV criterion value
  optimal_bandwidth <- bandwidth_range[which.min(lscv_values)]
  
  # Return the optimal bandwidth
  return(optimal_bandwidth)
}

mise_kde_old <- function(alltimes, status, h_vals, kernel = "triangular") {
  num <- length(alltimes)
  G_hat <- survfit(Surv(alltimes, 1 - status) ~ 1) # KM estimator for censoring survival
  scalier_fun <- stepfun(G_hat$time, c(1, G_hat$surv))
  
  # Function to compute leave-one-out MISE for a specific bandwidth
  leave_one_out_mise <- function(h) {
    squared_errors <- numeric(num) # To store the squared errors for each observation
    
    for (i in 1:num) {
      # Leave out the i-th observation
      alltimes_loo <- alltimes[-i]
      status_loo <- status[-i]
      uncensored_times <- alltimes_loo[status_loo == 1]
      
      # KDE estimation at the i-th observation, leaving it out
      if (kernel == "triangular") {
        kde_i <- sum(K_triang((uncensored_times - alltimes[i]) / h) / scalier_fun(uncensored_times)) / ((num - 1) * h)
      } else if (kernel == "gaussian") {
        kde_i <- sum(K_gaussian((uncensored_times - alltimes[i]) / h) / scalier_fun(uncensored_times)) / ((num - 1) * h)
      }
      
      # Compute the squared error for the i-th observation
      squared_errors[i] <- (kde_i - sum(kde_i) / ((num - 1) * h))^2 
    }
    
    # Mean of squared errors (MISE approximation)
    return(mean(squared_errors))
  }
  mise_errors <- sapply(h_vals, leave_one_out_mise)
  
  # Find the best bandwidth (the one with the lowest MISE)
  optimal_h <- h_vals[which.min(mise_errors)]
  return(list(optimal_h = optimal_h, mise_errors = mise_errors))
}

############################################################################################
# Other test of equality of quantiles
# Tang Jeong (BC test with contingency tables)

pval_Tang_Jeong = function(data1, Time1, data2, Time2, p, alpha = 0.05, return_stat = FALSE){
  n1 = length(data1$event)
  n2 = length(data2$event)
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  data2_PLOT = data.frame(data2)
  data2_PLOT$groups = rep(2, n2)
  
  DATA_PLOT = rbind(data1_PLOT, data2_PLOT)
  group_var = DATA_PLOT$groups
  
  result_1 = survQuantile(Time1, data1$event)
  result_2 = survQuantile(Time2, data2$event)
  
  #group_medians <- tapply(DATA_PLOT$x, DATA_PLOT$groups, median) 
  group_quant = tapply(DATA_PLOT$x, DATA_PLOT$groups, function(x) quantile(x, probs = c(p)))
  
  #fit = survfit(Surv(DATA_PLOT$x, 1-DATA_PLOT$event)~1)
  fit = survfit(Surv(DATA_PLOT$x, DATA_PLOT$event)~1)
  #overall_median = unname(median(fit))
  overall_quant = unname(quantile(fit, probs = p)$quantile)
  
  # Rank data
  DATA_PLOT$rank <- rank(DATA_PLOT$x)
  
  # Compute sum of ranks for each group
  rank_sums <- tapply(DATA_PLOT$rank, DATA_PLOT$groups, sum)
  
  # Calculate expected ranks and variances for each group
  n_groups <- table(DATA_PLOT$groups) # sample sizes per group
  n_times <- length(DATA_PLOT$x) # total number of times observed (=n=n1+n2)
  
  # Kaplan Meier survival curves
  F_1 = survfit(Surv(data1$x, data1$event)~1)
  #F_1 = survfit(Surv(data1$x, 1- data1$event)~1)
  scalier_fun_1 <- stepfun(F_1$time,c(1,F_1$surv))
  
  F_2 = survfit(Surv(data2$x, data2$event)~1)
  #F_2 = survfit(Surv(data2$x, 1- data2$event)~1)
  scalier_fun_2 <- stepfun(F_2$time,c(1,F_2$surv))
  #########
  # Tang and Jeong: contingency table
  #group 1
  n_11 = 0
  for(ind in 1:n1){
    if(data1$event[ind] == 1){
      if(data1$x[ind] > overall_quant){ # CHECA O TEMPO CENSURADO, NAO O TEMPO VERDADE
        n_11 = n_11 + 1
      }
      if(data1$x[ind] <= overall_quant){
        n_11 = n_11 
      }
    }
    if(data1$event[ind] == 0){
      if(data1$x[ind] >= overall_quant){
        n_11 = n_11 + 1
      }
      if(data1$x[ind] < overall_quant){
        n_11 = n_11+ scalier_fun_1(overall_quant)/scalier_fun_1( data1$x[ind])
      }
      
    }
  }
  
  n_12 = n1 - n_11
  
  # group 2
  n_21 = 0
  for(ind in 1:n2){
    if(data2$event[ind] == 1){
      if(data2$x[ind] > overall_quant){
        n_21 = n_21 + 1
      }
      if(data2$x[ind] <= overall_quant){
        n_21 = n_21 
      }
    }
    if(data2$event[ind] == 0){
      if(data2$x[ind] >= overall_quant){
        n_21 = n_21 + 1
      }
      if(data2$x[ind] < overall_quant){
        n_21 = n_21 + scalier_fun_2(overall_quant)/scalier_fun_2( data2$x[ind])
      }
      
    }
  }
  
  n_22 = n2 - n_21
  
  cont_table1 = data.frame(bigger = c(floor(n_11), floor(n_21)) , smaller = c(n1 - floor(n_11), n2 - floor(n_21)), total = c(n1, n2))
  cont_table2 = data.frame(bigger = c(floor(n_11), floor(n_21) + 1) , smaller = c(n1 - floor(n_11), n2 - (floor(n_21) + 1)), total = c(n1, n2)) 
  cont_table3 = data.frame(bigger = c(floor(n_11) + 1, floor(n_21) + 1) , smaller = c(n1 - (floor(n_11) + 1), n2 - (floor(n_21) + 1)), total = c(n1, n2)) 
  cont_table4 = data.frame(bigger = c(floor(n_11)+1, floor(n_21) ) , smaller = c(n1 - (floor(n_11)+1), n2 - floor(n_21) ), total = c(n1, n2)) 
  
  N = n1 + n2
  
  tables = list(cont_table1, cont_table2, cont_table3, cont_table4)
  Vs = c()
  for(tab in tables){
    mu_11 = (n1*sum(tab$bigger) )/N
    mu_12 = (n1*sum(tab$smaller) )/N
    mu_21 = (n2*sum(tab$bigger) )/N
    mu_22 = (n2*sum(tab$smaller) )/N
    V_i = (tab$bigger[1] -  mu_11)^2/mu_11 + (tab$smaller[1] -  mu_12)^2/mu_12 + (tab$bigger[2] -  mu_21)^2/mu_21 + (tab$bigger[2] -  mu_22)^2/mu_22
    Vs = c(Vs, V_i)
  }
  
  #U = sum(Vs) # follows a chi squared with 1 degree of freedom
  U = sum(0.25*Vs) # considering the weights
  
  critical_value = qchisq(1-alpha, 1)
  #critical_value = qchisq(alpha, 1)
  
  # Calculate the p-value
  p_value_TJ <- 1 - pchisq(U, 1)
  if(return_stat){
    return(list("pval" = p_value_TJ, "stat" = U))
  }
  return(p_value_TJ)
}


# KDE test

SurvQtest_kde <- function(data1, Time1, data2, Time2, p, h, ker = "triangular", alpha = 0.05) {
  # INPUT:
  # data should be sorted!!!!!
  # data1 (Required) = Random sample from F. It consists of two columns. First column contains the observed times. Second column contains censoring value( 1 = event, 0 = censoring ).
  # data2 (Required) = Random sample from G. It consists of two columns. First column contains the observed times. Second column contains censoring value( 1 = event, 0 = censoring ).
  # Time1 (Required) = Variable containing time to event or last follow-up in any units in data1.
  # Time2 (Required) = Variable containing time to event or last follow-up in any units in data2.
  # p (Required) = The probability of which the quantile would be tested. The null hypothesis is that the pth quantile of the F distribution is same as the pth quantile of the G distribution. .
  # h (Required) = Bandwidth adjustment term for ker density estimation.
  # ker (Default is triangular) = ker function for the density estimation. Can be either "triangular" or "gaussian"
  # alpha (Default is 0.05) = Type I error.
  # delta = fixed difference of quantiles under H1
  # OUTPUT:
  # result = Dataframe with NUM_F, NUM_G, XI_F, XI_G, DIFF, STD_ERR, LOWER, UPPER, CHISQ, PVALUE
  
  ##########
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  # In R, the equivalent functionality to the proc lifetest in SAS is provided by the survfit function from the survival package (compute KM)
  
  # Run lifetest on data1
  out1 <- survfit(Surv(x, event == 0) ~ x, data = data1)
  #out1 <- survfit(Surv(x, event == 0) ~ 1, data = data1)
  
  #out1_df <- data.frame(x = out1$time, event = data1$event, survival = out1$surv)
  
  
  # Run lifetest on data2
  out2 <- survfit(Surv(x, event == 0) ~ x, data = data2)
  #out2 <- survfit(Surv(x, event == 0) ~ 1, data = data2)
  #out2_df <- data.frame(x = out2$time, event = data2$event, survival = out2$surv)
  
  
  # Filter out1_df and out2_df
  #out1_df <- out1_df %>%
  #filter(x > 0 & event== 0)
  
  #out2_df <- out2_df %>%
  #filter(x > 0 & event == 0)
  
  # Calculate the minimum survival times for data1 and data2
  mymin3 <- min(out1$time)
  mymin4 <- min(out2$time)
  
  # Calculate m = 1 - max(mymin3, mymin4)
  m <- 1 - max(mymin3, mymin4)
  
  # Check if m < P and display an error message if true
  if (m < p) {
    cat('ERROR: Probability <P> is not appropriate. Try another one.\n')
    #stop("ERROR - detected in the input data to the function <SurvQtest>.")
  }
  
  # Sort data sets
  out1 <- arrange(data1, Time1)
  out2 <- arrange(data2, Time2)
  
  # sorting the data
  Tsort1 = sort(Time1, index.return = TRUE)
  alltimes1 = Tsort1$x
  status1 = data1$event[Tsort1$ix]
  
  Tsort2 = sort(Time2, index.return = TRUE)
  alltimes2 = Tsort2$x
  status2 = data2$event[Tsort2$ix]
  
  # Read data sets into matrices
  x <- as.matrix(out1)
  y <- as.matrix(out2)
  
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  xtmp1 <- matrix(0, n1, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n1)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n1)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n1) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  
  xtmp4 <- numeric(n1)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n1) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n1) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  # single numeric value representing the minimum index in the xtmp4 vector where the value is greater than or equal to the specified threshold p
  index1 <- min(which(xtmp4 >= p))
  
  # value from the first column of the matrix x at the row index1
  xxi <- x[index1, 1]
  
  # single numeric value representing the element from vector xtmp3 at position index1
  xhats <- xtmp3[index1]
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  xden <- (xtmp1[, 2] - xtmp1[, 1]) * xtmp1[, 2]
  
  xtmp6 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the division of the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n1) {
    if (xden[i] != 0) {
      xtmp6[i] <- xtmp1[i, 1] / xden[i]
    }
  }
  
  # Calculate phi: single numeric value
  phi <- (xhats^2) * sum((x[, 1] <= xxi) * xtmp6)
  
  
  # now we compute the equivalent vectors and matrix considering out 2 and y (instead of out1 and x)
  ytmp1 <- matrix(0, n2, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n2) {
    ytmp1[i, 1] <- as.numeric(y[i, 2] == 1)
    ytmp1[i, 2] <- sum(y[, 1] >= y[i, 1])
  }
  
  ytmp2 <- numeric(n2)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n2) {
    ytmp2[i] <- 1 - ytmp1[i, 1] / ytmp1[i, 2]
  }
  
  ytmp3 <- numeric(n2)
  for (i in 1:n2) {
    ytmp3[i] <- 1
    for (k in 1:i) {
      ytmp3[i] <- ytmp3[i] * ytmp2[k]
    }
  }
  
  ytmp4 <- numeric(n2)
  for (i in 1:n2) {
    ytmp4[i] <- 1 - ytmp3[i]
  }
  
  ytmp5 <- numeric(n2)
  ytmp5[1] <- 1 - ytmp4[1]
  for (i in 2:n2) {
    ytmp5[i] <- ytmp4[i] - ytmp4[i - 1]
  }
  
  # single numeric value representing the minimum index in the xtmp4 vector where the value is greater than or equal to the specified threshold p
  index2 <- min(which(ytmp4 >= p))
  
  # value from the first column of the matrix x at the row index1
  yxi <- y[index2, 1]
  
  # single numeric value representing the element from vector xtmp3 at position index1
  yhats <- ytmp3[index2]
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  yden <- (ytmp1[, 2] - ytmp1[, 1]) * ytmp1[, 2]
  
  ytmp6 <- numeric(n2)
  # vector of length n1
  # each element corresponds to the division of the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n2) {
    if (yden[i] != 0) {
      ytmp6[i] <- ytmp1[i, 1] / yden[i]
    }
  }
  
  
  # Calculate gamma: single numeric value
  gamma <- (yhats^2) * sum((y[, 1] <= yxi) * ytmp6)
  
  #####################################################################
  # Computations:
  # ker density estimator
  ########################
  #f = our_kde(xxi, alltimes1, status1, h, ker) 
  #g = our_kde(yxi, alltimes2 , status2, h, ker) 
  
  # with R's density estimation by kde
  kde_f = density(alltimes1, kernel = ker,adjust=1, bw = h)
  f = with(kde_f, approxfun(x, y)(xxi))
  kde_g = density(alltimes2, kernel = ker,adjust=1, bw = h)
  g = with(kde_g, approxfun(x, y)(yxi))
  
  phi = phi*n1
  gamma = gamma*n2
  
  delta <- xxi - yxi
  
  var_H0 <- phi /(mu*(f^2) ) + gamma /((1-mu)*(g^2))
  
  #stat <- delta^2/var_H0 # its squared, so compare to chisq
  stat = delta/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  # print result:
  #alpha <- ALPHA
  NUM_F <- n1
  NUM_G <- n2
  XI_F <- xxi
  XI_G <- yxi
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  #PVALUE <- 1 - pchisq(stat, df = 1) # because the stat is squared, we use chisq instead of gaussian
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = f, 'delta' = delta))
  
  result <- data.frame(
    NUM_F = NUM_F,
    NUM_G = NUM_G,
    XI_F = XI_F,
    XI_G = XI_G,
    delta = delta,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    STAT = stat, 
    PVALUE = PVALUE, 
    f = f,
    g = g,
    PHI = phi,
    GAMMA = gamma,
    POWER = POWER
  )
  
  return(result)
  #return(PVALUE)
}

###
# lin test

SurvQtestLin <- function(data1, Time1, data2, Time2, p, sd = 1, B=100, alpha = 0.05) {
  # INPUT:
  # data1 (Required) = Random sample from F. It consists of two columns. First column contains the observed times. Second column contains censoring value( 1 = event, 0 = censoring ).
  # data2 (Required) = Random sample from G. It consists of two columns. First column contains the observed times. Second column contains censoring value( 1 = event, 0 = censoring ).
  # Time1 (Required) = Variable containing time to event or last follow-up in any units in DATA1.
  # Time2 (Required) = Variable containing time to event or last follow-up in any units in DATA2.
  # p (Required) = The probability of which the quantile would be tested. The null hypothesis is that the pth quantile of the F distribution is same as the pth quantile of the G distribution. .
  # h (Required) = Bandwidth adjustment term for kernel density estimation.
  # kernel (Default is triangular) = Kernel function for the density estimation. Can be either "triangular" or "gaussian"
  # alpha (Default is 0.05) = Type I error.
  # delta = fixed difference of quantiles under H1
  # OUTPUT:
  # result = Dataframe with NUM_F, NUM_G, XI_F, XI_G, DIFF, STD_ERR, LOWER, UPPER, CHISQ, PVALUE
  
  ##########
  # In R, the equivalent functionality to the proc lifetest in SAS is provided by the survfit function from the survival package (compute KM)
  
  # Run lifetest on data1
  out1 <- survfit(Surv(x, event == 0) ~ x, data = data1)
  #out1 <- survfit(Surv(x, event == 0) ~ 1, data = data1)
  
  #out1_df <- data.frame(x = out1$time, event = data1$event, survival = out1$surv)
  
  
  # Run lifetest on data2
  out2 <- survfit(Surv(x, event == 0) ~ x, data = data2)
  #out2 <- survfit(Surv(x, event == 0) ~ 1, data = data2)
  #out2_df <- data.frame(x = out2$time, event = data2$event, survival = out2$surv)
  
  
  # Filter out1_df and out2_df
  #out1_df <- out1_df %>%
  #filter(x > 0 & event== 0)
  
  #out2_df <- out2_df %>%
  #filter(x > 0 & event == 0)
  
  # Calculate the minimum survival times for data1 and data2
  mymin3 <- min(out1$time)
  mymin4 <- min(out2$time)
  
  # Calculate m = 1 - max(mymin3, mymin4)
  m <- 1 - max(mymin3, mymin4)
  
  # Check if m < P and display an error message if true
  if (m < p) {
    cat('ERROR: Probability <P> is not appropriate. Try another one.\n')
    #stop("ERROR - detected in the input data to the function <SurvQtest>.")
  }
  
  # Sort data sets
  out1 <- arrange(data1, Time1)
  out2 <- arrange(data2, Time2)
  
  # sorting the data
  Tsort1 = sort(Time1, index.return = TRUE)
  alltimes1 = Tsort1$x
  status1 = data1$event[Tsort1$ix]
  
  Tsort2 = sort(Time2, index.return = TRUE)
  alltimes2 = Tsort2$x
  status2 = data2$event[Tsort2$ix]
  
  # Read data sets into matrices
  x <- as.matrix(out1)
  y <- as.matrix(out2)
  
  n1 <- nrow(x)
  n2 <- nrow(y)
  
  xtmp1 <- matrix(0, n1, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n1)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n1)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n1) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  
  xtmp4 <- numeric(n1)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n1) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n1) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  # single numeric value representing the minimum index in the xtmp4 vector where the value is greater than or equal to the specified threshold p
  index1 <- min(which(xtmp4 >= p))
  
  # value from the first column of the matrix x at the row index1
  xxi <- x[index1, 1]
  
  # single numeric value representing the element from vector xtmp3 at position index1
  xhats <- xtmp3[index1]
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  xden <- (xtmp1[, 2] - xtmp1[, 1]) * xtmp1[, 2]
  
  xtmp6 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the division of the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n1) {
    if (xden[i] != 0) {
      xtmp6[i] <- xtmp1[i, 1] / xden[i]
    }
  }
  
  # Calculate phi: single numeric value
  phi <- (xhats^2) * sum((x[, 1] <= xxi) * xtmp6)
  
  # now we compute the equivalent vectors and matrix considering out 2 and y (instead of out1 and x)
  
  ytmp1 <- matrix(0, n2, 2)
  # matrix of size n2 by 2
  # similar to xtmp1 but for out2 and y (instead of out1 and x)
  for (i in 1:n2) {
    ytmp1[i, 1] <- as.numeric(y[i, 2] == 1)
    ytmp1[i, 2] <- sum(y[, 1] >= y[i, 1])
  }
  
  ytmp2 <- numeric(n2)
  for (i in 1:n2) {
    ytmp2[i] <- 1 - ytmp1[i, 1] / ytmp1[i, 2]
  }
  
  ytmp3 <- numeric(n2)
  for (i in 1:n2) {
    ytmp3[i] <- 1
    for (k in 1:i) {
      ytmp3[i] <- ytmp3[i] * ytmp2[k]
    }
  }
  
  ytmp4 <- numeric(n2)
  for (i in 1:n2) {
    ytmp4[i] <- 1 - ytmp3[i]
  }
  
  ytmp5 <- numeric(n2)
  ytmp5[1] <- 1 - ytmp4[1]
  for (i in 2:n2) {
    ytmp5[i] <- ytmp4[i] - ytmp4[i - 1]
  }
  
  index2 <- min(which(ytmp4 >= p))
  
  yxi <- y[index2, 1]
  yhats <- ytmp3[index2]
  
  yden <- (ytmp1[, 2] - ytmp1[, 1]) * ytmp1[, 2]
  
  ytmp6 <- numeric(n2)
  for (i in 1:n2) {
    if (yden[i] != 0) {
      ytmp6[i] <- ytmp1[i, 1] / yden[i]
    }
  }
  
  # Calculate gamma: single numeric value
  gamma <- (yhats^2) * sum((y[, 1] <= yxi) * ytmp6)
  
  #single numeric value representing the index in the vector where the value is greater than or equal to the specified threshold p
  index3 <- min(which(xtmp4 >= p)) # antes era q
  
  #retrieves the value from the x vector at the index specified by index3, in the first column.
  xq_xi <- x[index3, 1]
  
  
  #####################################################################
  # Computations:
  # Lin density estimator
  ########################
  # Kaplan Meier estimator of the survival
  F_hat = survfit(Surv(x, event == 1) ~ 1, data = data1) # KM estimator of the survival 
  F_scalier_fun <- stepfun(F_hat$time,c(1,F_hat$surv))
  
  # distribution function
  F_dist = function(x){
    1- F_scalier_fun(x) 
  }
  
  G_hat = survfit(Surv(x, event == 1) ~ 1, data = data2) # KM estimator of the survival 
  G_scalier_fun <- stepfun(G_hat$time,c(1,G_hat$surv))
  
  # distribution function
  G_dist = function(x){
    1- G_scalier_fun(x) 
  }
  
  ########################
  # density estimator
  Zi_vec1 = rnorm(n = B, mean = 0, sd = sd)
  Zi_vec2 = rnorm(n = B, mean = 0, sd = sd)
  Yi_vec1 = sqrt(n1)*( F_dist( xxi + Zi_vec1/sqrt(n1) ) - p )
  #Yi_vec2 = sqrt(n2)*( G_dist( yxi + Zi_vec2/sqrt(n2) ) - p )
  Yi_vec2 = sqrt(n2)*( G_dist( yxi + Zi_vec1/sqrt(n2) ) - p )
  
  
  # Zi_vec1 = c()
  # Yi_vec1 = c()
  # 
  # Zi_vec2 = c()
  # Yi_vec2 = c()
  # for(it in 1:B){
  #   Zi1 = rnorm(n = 1, mean = 0, sd = sd)
  #   Yi1 = sqrt(n1)*( F_dist( xxi + Zi1/sqrt(n1) ) - p )
  #   
  #   Zi2 = rnorm(n = 1, mean = 0, sd = sd)
  #   Yi2 = sqrt(n2)*( G_dist( yxi + Zi2/sqrt(n2) ) - p )
  #   
  #   Zi_vec1[it] = Zi1
  #   Yi_vec1[it] = Yi1
  #   
  #   Zi_vec2[it] = Zi2
  #   Yi_vec2[it] = Yi2
  # }
  # 
  
  model1 = lm( Yi_vec1 ~ (Zi_vec1) -1 )
  model2 = lm( Yi_vec2 ~ (Zi_vec1) -1 )
  #model2 = lm( Yi_vec2 ~ (Zi_vec2) -1 )
  f = unname(coef(model1))
  g = unname(coef(model2))
  
  #model1 = gam(Yi_vec1 ~ Zi_vec1)
  #model2 = gam(Yi_vec2 ~ Zi_vec1)
  #f = -unname(coef(model1))[1]
  #g = -unname(coef(model2))[1]
  
  phi = phi*n1
  gamma = gamma*n2
  n = n1+n2 
  mu = n1/n
  
  delta <- (xxi - yxi)
  psi <- phi /(mu *(f * f)) + gamma /((1-mu)* (g * g))
  #stat <- (delta^2) / psi # its squared, so compare to chisq
  
  stat = delta/sqrt(psi)
  
  stat = sqrt(n)*stat
  
  
  # print result:
  #alpha <- ALPHA
  NUM_F <- n1
  NUM_G <- n2
  XI_F <- xxi
  XI_G <- yxi
  DIFF <- xxi - yxi
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(psi)
  LOWER <- DIFF - z * sqrt(psi)
  UPPER <- DIFF + z * sqrt(psi)
  STAT <- stat
  #PVALUE <- 1 - pchisq(stat, df = 1) # because the stat is squared, we use chisq instead of gaussian
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  POWER <- 1-pnorm(z - sqrt((n1+n2)/psi)*delta ) + pnorm(-z - sqrt((n1+n2)/psi)*delta )
  
  result <- data.frame(
    NUM_F = NUM_F,
    NUM_G = NUM_G,
    XI_F = XI_F,
    XI_G = XI_G,
    DIFF = DIFF,
    STD_ERR = STD_ERR[1], # multiple values
    LOWER = LOWER[1], # multiple values
    UPPER = UPPER[1], # multiple values
    STAT = STAT[1], # multiple values
    PVALUE = PVALUE[1], # multiple values
    f = f,
    g=g,
    PHI = phi,
    GAMMA = gamma,
    POWER = POWER
  )
  
  return(result)
  #return(PVALUE)
}


###################
# Multivariate test
######################
SurvQtestLin_Multi <- function(data1, Time1, data2, Time2, p1, p2, sd = 1, B=100, alpha = 0.05) {
  
  alltimes1_unsort = data1$x
  status1_unsort = data1$event
  alltimes2_unsort = data2$x
  status2_unsort = data2$event
  
  km_fit_F <- survfit(Surv(alltimes1_unsort, status1_unsort) ~ 1)
  quantiles_F <- quantile(km_fit_F, probs = c(p1, p2))$quantile
  
  km_fit_G <- survfit(Surv(alltimes2_unsort, status2_unsort) ~ 1)
  quantiles_G <- quantile(km_fit_G, probs = c(p1, p2))$quantile
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  n = n1+n2
  
  mu = n1/n
  
  min_quant_F = min(unname(quantiles_F))
  min_quant_G = min(unname(quantiles_G))
  
  x = arrange(data1, Time1)
  y = arrange(data2, Time2)
  
  xtmp1 <- matrix(0, n1, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n1)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n1)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n1) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  
  xtmp4 <- numeric(n1)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n1) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n1) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  xden <- (xtmp1[, 2] - xtmp1[, 1]) * xtmp1[, 2]
  
  xtmp6 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the division of the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n1) {
    if (xden[i] != 0) {
      xtmp6[i] <- xtmp1[i, 1] / xden[i]
    }
  }
  
  
  Sigma_F_11 = ((1-p1)^2) * sum((x[, 1] <= quantiles_F[[1]]) * xtmp6)
  Sigma_F_22 = ((1-p2)^2) * sum((x[, 1] <= quantiles_F[[2]]) * xtmp6)
  Sigma_F_12 = (1-p1)*(1-p2) * sum((x[, 1] <=  min_quant_F) * xtmp6)
  Sigma_F_21 = Sigma_F_12
  
  
  ytmp1 <- matrix(0, n2, 2)
  # matrix of size n2 by 2
  # similar to xtmp1 but for out2 and y (instead of out1 and x)
  for (i in 1:n2) {
    ytmp1[i, 1] <- as.numeric(y[i, 2] == 1)
    ytmp1[i, 2] <- sum(y[, 1] >= y[i, 1])
  }
  
  ytmp2 <- numeric(n2)
  for (i in 1:n2) {
    ytmp2[i] <- 1 - ytmp1[i, 1] / ytmp1[i, 2]
  }
  
  ytmp3 <- numeric(n2)
  for (i in 1:n2) {
    ytmp3[i] <- 1
    for (k in 1:i) {
      ytmp3[i] <- ytmp3[i] * ytmp2[k]
    }
  }
  
  ytmp4 <- numeric(n2)
  for (i in 1:n2) {
    ytmp4[i] <- 1 - ytmp3[i]
  }
  
  ytmp5 <- numeric(n2)
  ytmp5[1] <- 1 - ytmp4[1]
  for (i in 2:n2) {
    ytmp5[i] <- ytmp4[i] - ytmp4[i - 1]
  }
  
  index2 <- min(which(ytmp4 >= p))
  
  yden <- (ytmp1[, 2] - ytmp1[, 1]) * ytmp1[, 2]
  
  ytmp6 <- numeric(n2)
  for (i in 1:n2) {
    if (yden[i] != 0) {
      ytmp6[i] <- ytmp1[i, 1] / yden[i]
    }
  }
  
  Sigma_G_11 = ((1-p1)^2) * sum((y[, 1] <= quantiles_G[[1]]) * ytmp6)
  Sigma_G_22 = ((1-p2)^2) * sum((y[, 1] <= quantiles_G[[2]]) * ytmp6)
  Sigma_G_12 = (1-p1)*(1-p2) * sum((y[, 1] <=  min_quant_G) * ytmp6)
  Sigma_G_21 = Sigma_G_12
  
  
  #####################################################################
  # Computations:
  # Lin density estimator
  ########################
  # Kaplan Meier estimator of the survival
  F_hat = survfit(Surv(x, event == 1) ~ 1, data = data1) # KM estimator of the survival 
  F_scalier_fun <- stepfun(F_hat$time,c(1,F_hat$surv))
  
  # distribution function
  F_dist = function(x){
    1- F_scalier_fun(x) 
  }
  
  G_hat = survfit(Surv(x, event == 1) ~ 1, data = data2) # KM estimator of the survival 
  G_scalier_fun <- stepfun(G_hat$time,c(1,G_hat$surv))
  
  # distribution function
  G_dist = function(x){
    1- G_scalier_fun(x) 
  }
  
  ########################
  # density estimator
  Zi_vec1 = rnorm(n = B, mean = 0, sd = sd)
  
  Yi_vec1_f = sqrt(n1)*( F_dist( quantiles_F[[1]] + Zi_vec1/sqrt(n1) ) - p1 )
  Yi_vec2_f = sqrt(n1)*( F_dist( quantiles_F[[2]] + Zi_vec1/sqrt(n1) ) - p2 )
  
  Yi_vec1_g = sqrt(n2)*( G_dist( quantiles_G[[1]] + Zi_vec1/sqrt(n2) ) - p1 )
  Yi_vec2_g = sqrt(n2)*( G_dist( quantiles_G[[2]] + Zi_vec1/sqrt(n2) ) - p2 )
  
  
  model1_f = lm( Yi_vec1_f ~ (Zi_vec1) -1 )
  model2_f = lm( Yi_vec2_f ~ (Zi_vec1) -1 )
  
  model1_g = lm( Yi_vec1_g ~ (Zi_vec1) -1 )
  model2_g = lm( Yi_vec2_g ~ (Zi_vec1) -1 )
  
  
  f_1 = unname(coef(model1_f))
  f_2 = unname(coef(model2_f))
  g_1 = unname(coef(model1_g))
  g_2 = unname(coef(model2_g))
  
  
  #model1 = gam(Yi_vec1 ~ Zi_vec1)
  #model2 = gam(Yi_vec2 ~ Zi_vec1)
  #f = -unname(coef(model1))[1]
  #g = -unname(coef(model2))[1]
  
  
  Sigma_G = matrix(
    c(Sigma_G_11/((1-mu)*g_1^2), Sigma_G_12/((1-mu)*g_1*g_2 ), 
      Sigma_G_21/((1-mu)*g_1*g_2 ), Sigma_G_22/((1-mu)*(g_2^2))), 
    nrow = 2, byrow = TRUE
  )
  
  Sigma_F = matrix(
    c(Sigma_F_11/(mu*f_1^2), Sigma_F_12/(mu*f_1*f_2), 
      Sigma_F_21/(mu*f_1*f_2), Sigma_F_22/(mu*f_2^2)), 
    nrow = 2, byrow = TRUE
  )
  
  upsilon = Sigma_F+Sigma_G
  upsilon_inv = solve(upsilon)
  
  Zn = matrix(
    c(quantiles_F[[1]] - quantiles_G[[1]],
      quantiles_F[[2]] - quantiles_G[[2]])
  )
  
  #Zn = sqrt(n)*Zn
  
  stat = t(Zn) %*% upsilon_inv %*% Zn ##### VER SE DE FATO INVERTE!!
  
  
  #stat = sqrt(n)*stat
  
  # print result:
  #alpha <- ALPHA
  NUM_F <- n1
  NUM_G <- n2
  
  STAT <- stat
  #PVALUE <- 1 - pchisq(stat, df = 1) # because the stat is squared, we use chisq instead of gaussian
  PVALUE <- 1 - pchisq(stat, 2)
  #POWER <- 1-pnorm(z - sqrt((n1+n2)/psi)*delta ) + pnorm(-z - sqrt((n1+n2)/psi)*delta )
  
  result <- data.frame(
    NUM_F = NUM_F,
    NUM_G = NUM_G,
    STAT = STAT[1], # multiple values
    PVALUE = PVALUE[1]
  )
  print(upsilon_inv)
  
  print(Zn)
  
  return(result)
  #return(PVALUE)
}


SurvQtestLin_Multi_kde <- function(data1, Time1, data2, Time2, p1, p2, ker = "triangular", h = 0.1, alpha = 0.05) {
  
  alltimes1_unsort = data1$x
  status1_unsort = data1$event
  alltimes2_unsort = data2$x
  status2_unsort = data2$event
  
  km_fit_F <- survfit(Surv(alltimes1_unsort, status1_unsort) ~ 1)
  quantiles_F <- quantile(km_fit_F, probs = c(p1, p2))$quantile
  
  km_fit_G <- survfit(Surv(alltimes2_unsort, status2_unsort) ~ 1)
  quantiles_G <- quantile(km_fit_G, probs = c(p1, p2))$quantile
  
  n1 <- nrow(data1)
  n2 <- nrow(data2)
  
  n = n1+n2
  
  mu = n1/n
  
  min_quant_F = min(unname(quantiles_F))
  min_quant_G = min(unname(quantiles_G))
  
  x = arrange(data1, Time1)
  y = arrange(data2, Time2)
  
  xtmp1 <- matrix(0, n1, 2) 
  # 2 column matrix: 
  #the 1st indicates whether the value in the 2nd column of the corresponding row in out1 is equal to 1
  #the 2nd represents the count of rows in out1 where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp1[i, 1] <- as.numeric(x[i, 2] == 1)
    xtmp1[i, 2] <- sum(x[, 1] >= x[i, 1])
  }
  
  xtmp2 <- numeric(n1)
  # vector of length n1 (=number of rows in out1), each element corresponds to a row in the out1 dataset
  # computes proportion of times the value in the second column of the corresponding row of out1 is not equal to 1, among the times where the value in the first column is greater than or equal to the value in the corresponding row
  for (i in 1:n1) {
    xtmp2[i] <- 1 - xtmp1[i, 1] / xtmp1[i, 2]
  }
  
  xtmp3 <- numeric(n1)
  # vector of length n1
  #vector where each element is the cumulative product of specific elements from the xtmp2 vector
  for (i in 1:n1) {
    xtmp3[i] <- 1
    for (k in 1:i) {
      xtmp3[i] <- xtmp3[i] * xtmp2[k]
    }
  }
  
  xtmp4 <- numeric(n1)
  # vector of length n1
  # each element has the difference between 1 and the corresponding element of xtmp3
  for (i in 1:n1) {
    xtmp4[i] <- 1 - xtmp3[i]
  }
  
  xtmp5 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the difference between consecutive elements from xtmp4
  xtmp5[1] <- 1 - xtmp4[1]
  for (i in 2:n1) {
    xtmp5[i] <- xtmp4[i] - xtmp4[i - 1]
  }
  
  
  # vector where each element is: (difference between 2nd and 1st columns of xtmp1 for each row) times (values from the 2nd column of xtmp1 for each row)
  xden <- (xtmp1[, 2] - xtmp1[, 1]) * xtmp1[, 2]
  
  xtmp6 <- numeric(n1)
  # vector of length n1
  # each element corresponds to the division of the element in the 1st column of xtmp1 by the corresponding element from xden
  for (i in 1:n1) {
    if (xden[i] != 0) {
      xtmp6[i] <- xtmp1[i, 1] / xden[i]
    }
  }
  
  
  Sigma_F_11 = ((1-p1)^2) * sum((x[, 1] <= quantiles_F[[1]]) * xtmp6)
  Sigma_F_22 = ((1-p2)^2) * sum((x[, 1] <= quantiles_F[[2]]) * xtmp6)
  Sigma_F_12 = (1-p1)*(1-p2) * sum((x[, 1] <=  min_quant_F) * xtmp6)
  Sigma_F_21 = Sigma_F_12
  
  
  ytmp1 <- matrix(0, n2, 2)
  # matrix of size n2 by 2
  # similar to xtmp1 but for out2 and y (instead of out1 and x)
  for (i in 1:n2) {
    ytmp1[i, 1] <- as.numeric(y[i, 2] == 1)
    ytmp1[i, 2] <- sum(y[, 1] >= y[i, 1])
  }
  
  ytmp2 <- numeric(n2)
  for (i in 1:n2) {
    ytmp2[i] <- 1 - ytmp1[i, 1] / ytmp1[i, 2]
  }
  
  ytmp3 <- numeric(n2)
  for (i in 1:n2) {
    ytmp3[i] <- 1
    for (k in 1:i) {
      ytmp3[i] <- ytmp3[i] * ytmp2[k]
    }
  }
  
  ytmp4 <- numeric(n2)
  for (i in 1:n2) {
    ytmp4[i] <- 1 - ytmp3[i]
  }
  
  ytmp5 <- numeric(n2)
  ytmp5[1] <- 1 - ytmp4[1]
  for (i in 2:n2) {
    ytmp5[i] <- ytmp4[i] - ytmp4[i - 1]
  }
  
  index2 <- min(which(ytmp4 >= p))
  
  yden <- (ytmp1[, 2] - ytmp1[, 1]) * ytmp1[, 2]
  
  ytmp6 <- numeric(n2)
  for (i in 1:n2) {
    if (yden[i] != 0) {
      ytmp6[i] <- ytmp1[i, 1] / yden[i]
    }
  }
  
  Sigma_G_11 = ((1-p1)^2) * sum((y[, 1] <= quantiles_G[[1]]) * ytmp6)
  Sigma_G_22 = ((1-p2)^2) * sum((y[, 1] <= quantiles_G[[2]]) * ytmp6)
  Sigma_G_12 = (1-p1)*(1-p2) * sum((y[, 1] <=  min_quant_G) * ytmp6)
  Sigma_G_21 = Sigma_G_12
  
  
  #####################################################################
  # Computations:
  # Lin density estimator
  ########################
  # Kaplan Meier estimator of the survival
  F_hat = survfit(Surv(x, event == 1) ~ 1, data = data1) # KM estimator of the survival 
  F_scalier_fun <- stepfun(F_hat$time,c(1,F_hat$surv))
  
  # distribution function
  F_dist = function(x){
    1- F_scalier_fun(x) 
  }
  
  G_hat = survfit(Surv(x, event == 1) ~ 1, data = data2) # KM estimator of the survival 
  G_scalier_fun <- stepfun(G_hat$time,c(1,G_hat$surv))
  
  # distribution function
  G_dist = function(x){
    1- G_scalier_fun(x) 
  }
  
  ########################
  # density estimator
  
  kde_f = density(x$x, kernel = ker,adjust=1, bw = h)
  f_1 = with(kde_f, approxfun(x, y)(quantiles_F[[1]]))
  f_2 = with(kde_f, approxfun(x, y)(quantiles_F[[2]]))
  kde_g = density(y$x, kernel = ker,adjust=1, bw = h)
  g_1 = with(kde_g, approxfun(x, y)(quantiles_G[[1]]))
  g_2 = with(kde_g, approxfun(x, y)(quantiles_G[[2]]))
  
  #model1 = gam(Yi_vec1 ~ Zi_vec1)
  #model2 = gam(Yi_vec2 ~ Zi_vec1)
  #f = -unname(coef(model1))[1]
  #g = -unname(coef(model2))[1]
  
  
  Sigma_G = matrix(
    c(Sigma_G_11/((1-mu)*g_1^2), Sigma_G_12/((1-mu)*g_1*g_2 ), 
      Sigma_G_21/((1-mu)*g_1*g_2 ), Sigma_G_22/((1-mu)*(g_2^2))), 
    nrow = 2, byrow = TRUE
  )
  
  Sigma_F = matrix(
    c(Sigma_F_11/(mu*f_1^2), Sigma_F_12/(mu*f_1*f_2), 
      Sigma_F_21/(mu*f_1*f_2), Sigma_F_22/(mu*f_2^2)), 
    nrow = 2, byrow = TRUE
  )
  
  upsilon = Sigma_F+Sigma_G
  upsilon_inv = solve(upsilon)
  
  Zn = matrix(
    c(quantiles_F[[1]] - quantiles_G[[1]],
      quantiles_F[[2]] - quantiles_G[[2]])
  )
  
  #Zn = sqrt(n)*Zn
  
  stat = t(Zn) %*% upsilon_inv %*% Zn ##### VER SE DE FATO INVERTE!!
  
  
  #stat = sqrt(n)*stat
  
  # print result:
  #alpha <- ALPHA
  NUM_F <- n1
  NUM_G <- n2
  
  STAT <- stat
  #PVALUE <- 1 - pchisq(stat, df = 1) # because the stat is squared, we use chisq instead of gaussian
  PVALUE <- 1 - pchisq(stat, 2)
  #POWER <- 1-pnorm(z - sqrt((n1+n2)/psi)*delta ) + pnorm(-z - sqrt((n1+n2)/psi)*delta )
  
  result <- data.frame(
    NUM_F = NUM_F,
    NUM_G = NUM_G,
    STAT = STAT[1], # multiple values
    PVALUE = PVALUE[1]
  )
  print(upsilon_inv)
  
  print(Zn)
  
  return(result)
  #return(PVALUE)
}

#################
# planning
SurvQtestExp_new <- function(n1, n2, k1, k2, lamb_group1, p, delta, cens_type = 'unif', cens_rate = 0.3, alpha = 0.05, return_data = FALSE) {
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  data1 <- data.frame("x" = alltimes1_unsort , "event" = status1_unsort) #truetimes0_unsort
  Time1 <- alltimes1_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  quantile_group1 = -log(1-p)/lamb_group1
  quantile_group2 = -delta + quantile_group1
  
  lamb_group2 = -log(1-p)/quantile_group2
  
  truetimes2_unsort <- rexp(n2, rate = lamb_group2)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  data2 <- data.frame("x" = alltimes2_unsort, "event" = status2_unsort)
  Time2 <- alltimes2_unsort
  
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  data2_PLOT = data.frame(data2)
  data2_PLOT$groups = rep(2, n2)
  
  DATA_PLOT = rbind(data1_PLOT, data2_PLOT)
  
  #cat("censoring rate:", 1-(sum(data1_PLOT$event)/ length(data1_PLOT$event)) )
  
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  
  fun_integrand = function(rate, x){
    #return(rate*x*exp(rate*x))
    return(rate*exp(rate*x))
  }
  
  # estimating phi and gamma numerically
  #phi = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group1, rate = lamb_group1)$val # they are the same as computing with the 'true_phi' function
  #gamma = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group2, rate = lamb_group2)$val
  
  # computing explicit phi and gamma
  if(cens_type == 'unif'){
    if(quantile_group1 <k2){
      phi = (1-p)^2*(exp(lamb_group1*quantile_group1) - 1) 
      gamma = (1-p)^2*(exp(lamb_group2*quantile_group2) - 1) 
    }
    else{
      cat('ERROR: F-1(p) should be less than k2 which is minimum of the uniform censoring')
      stop("ERROR")
    }
  }
  
  else if(cens_type == 'exp'){
    phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
    gamma = (1-p)^2*(lamb_group2/(lamb_group2+cens_rate))*(exp((lamb_group2+cens_rate)*quantile_group2) -1)
  }
  
  
  var_H0 = phi/(mu*(lamb_group1 * exp(-lamb_group1 * quantile_group1))^2)  + gamma/((1-mu)* (lamb_group2 * exp(-lamb_group2 * quantile_group2))^2)
  
  stat = (delta)/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  #PVALUE <- 1 - pchisq(stat^2, df = 1) # because the stat is squared, we use chisq instead of gaussian
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = lamb_group1 * exp(-lamb_group1 * quantile_group1), 'delta' = delta))
  
  result <- data.frame(
    NUM_F = n1,
    NUM_G = n2,
    XI_F = quantile_group1,
    XI_G = quantile_group2,
    delta = delta,
    p = p,
    STAT = stat,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    PVALUE = PVALUE, 
    POWER = POWER
  )
  
  if(return_data){
    return(list(result, DATA_PLOT))
  }
  return(result)
  
}


SurvQtestExp_empirical <- function(n1, n2, k1, k2, lamb_group1, lamb_group2, p, cens_type = 'unif', cens_rate = 0.3, alpha = 0.05, return_data = FALSE) {
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  data1 <- data.frame("x" = alltimes1_unsort , "event" = status1_unsort) #truetimes0_unsort
  Time1 <- alltimes1_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  quantile_group1 = find_quantile(data1)
  
  #quantile_group1 = -log(1-p)/lamb_group1
  #quantile_group2 = -delta + quantile_group1
  
  #lamb_group2 = -log(1-p)/quantile_group2
  
  truetimes2_unsort <- rexp(n2, rate = lamb_group2)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  data2 <- data.frame("x" = alltimes2_unsort, "event" = status2_unsort)
  Time2 <- alltimes2_unsort
  
  quantile_group2 = find_quantile(data2)
  delta = quantile_group1 - quantile_group2
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  data2_PLOT = data.frame(data2)
  data2_PLOT$groups = rep(2, n2)
  
  DATA_PLOT = rbind(data1_PLOT, data2_PLOT)
  
  #cat("censoring rate:", 1-(sum(data1_PLOT$event)/ length(data1_PLOT$event)) )
  
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  
  fun_integrand = function(rate, x){
    #return(rate*x*exp(rate*x))
    return(rate*exp(rate*x))
  }
  
  # estimating phi and gamma numerically
  #phi = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group1, rate = lamb_group1)$val # they are the same as computing with the 'true_phi' function
  #gamma = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group2, rate = lamb_group2)$val
  
  # computing explicit phi and gamma
  if(cens_type == 'unif'){
    if(quantile_group1 <k2){
      phi = (1-p)^2*(exp(lamb_group1*quantile_group1) - 1) 
      gamma = (1-p)^2*(exp(lamb_group2*quantile_group2) - 1) 
    }
    else{
      cat('ERROR: F-1(p) should be less than k2 which is minimum of the uniform censoring')
      stop("ERROR")
    }
  }
  
  else if(cens_type == 'exp'){
    phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
    gamma = (1-p)^2*(lamb_group2/(lamb_group2+cens_rate))*(exp((lamb_group2+cens_rate)*quantile_group2) -1)
  }
  
  
  var_H0 = phi/(mu*(lamb_group1 * exp(-lamb_group1 * quantile_group1))^2)  + gamma/((1-mu)* (lamb_group2 * exp(-lamb_group2 * quantile_group2))^2)
  
  stat = (delta)/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  #PVALUE <- 1 - pchisq(stat^2, df = 1) # because the stat is squared, we use chisq instead of gaussian
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = lamb_group1 * exp(-lamb_group1 * quantile_group1), 'delta' = delta))
  
  result <- data.frame(
    NUM_F = n1,
    NUM_G = n2,
    XI_F = quantile_group1,
    XI_G = quantile_group2,
    delta = delta,
    p = p,
    STAT = stat,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    PVALUE = PVALUE, 
    POWER = POWER
  )
  
  if(return_data){
    return(list(result, DATA_PLOT))
  }
  return(result)
  
}

SurvQtestPieceExp_new <- function(n1, n2, k1, k2, lamb_group1, tcut, p, delta, cens_type = 'unif', cens_rate = 0.3, alpha = 0.05, late = TRUE) {
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  Tsort1 = sort(alltimes1_unsort, index.return = TRUE)
  alltimes1 = Tsort1$x
  status1 = alltimes1[Tsort1$ix]
  
  data1 <- data.frame("x" = alltimes1, "event" = status1) #truetimes0_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  quantile_group1 = -log(1-p)/lamb_group1
  quantile_group2 = -delta + quantile_group1
  
  if(late){
    lamb_begin = lamb_group1
    #if(quantile_group2 < tcut){
    if(p>=1-exp(-lamb_begin*tcut) & tcut-quantile_group2!=0){
      lamb_group2 = (log(1-p)+lamb_begin*tcut)/(tcut-quantile_group2)
      if(lamb_group2<0){
        lamb_group2 = -log(1-p)/quantile_group2
      }
    }
    else{
      print('error')
      lamb_group2 = -log(1-p)/quantile_group2
      if(lamb_group2<0){
        lamb_group2 = (log(1-p)+lamb_begin*tcut)/(tcut-quantile_group2)
      }
    }
    lamb_end = lamb_group2
  }
  else{
    lamb_end = lamb_group1
    lamb_begin_candidate1 = -log(1-p)/quantile_group2
    lamb_begin_candidate2 = (lamb_end*(tcut - quantile_group2) - log(1-p))/tcut
    if(p<1-exp(-lamb_begin_candidate1*tcut)){
      lamb_group2 = lamb_begin_candidate1
    }
    else{
      lamb_group2 = lamb_begin_candidate2
    }
    lamb_begin = lamb_group2
  }
  
  rate = c(lamb_begin, lamb_end)
  truetimes2_unsort = rsurv(n2, cuts = c(tcut), alpha = rate)
  
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n2, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n2, rate = cens_rate)
  }
  
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  Tsort2 = sort(alltimes2_unsort, index.return = TRUE)
  alltimes2 = Tsort2$x
  status2 = alltimes2[Tsort2$ix]
  
  data2 <- data.frame("x" = alltimes2, "event" = status2) #truetimes0_unsort
  
  Time2 <- alltimes2_unsort
  
  
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  
  #fun_integrand = function(rate, x){
  #return(rate*x*exp(rate*x))
  #}
  
  # estimating phi and gamma numerically
  #phi = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group1, rate = lamb_group1)$val # they are the same as computing with the 'true_phi' function
  #gamma = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group2, rate = lamb_group2)$val
  
  # computing explicit phi and gamma
  if(cens_type == 'unif'){
    if(quantile_group1 <k2){
      phi = (1-p)^2*(exp(lamb_group1*quantile_group1) - 1)
      cat('ERROR: Add computation for gamma')
      stop("ERROR")
    }
  }
  
  else if(cens_type == 'exp'){
    phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
    #gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) + (1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin-lamb_end)*tcut))*(exp((lamb_end+cens_rate)*(quantile_group2 - tcut)))
    gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) + (1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin+cens_rate)*tcut))*(exp((lamb_end+cens_rate)*(quantile_group2 - tcut)))
    
  }
  
  
  f = lamb_group1*exp(-lamb_group1 * quantile_group1)
  g = g_piecewise(tcut, lamb_begin, lamb_end, quantile_group2)
  
  var_H0 = phi/(mu*(f^2))  + gamma/((1-mu)* g^2)
  
  stat = (delta)/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  #PVALUE <- 1 - pchisq(stat^2, df = 1) # because the stat is squared, we use chisq instead of gaussian
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = f, 'g' = g, 'gamma' = gamma, 'delta' = delta))
  
  result <- data.frame(
    NUM_F = n1,
    NUM_G = n2,
    XI_F = quantile_group1,
    XI_G = quantile_group2,
    delta = delta,
    p = p,
    LAMB_GROUP2 = lamb_group2,
    STAT = stat,
    g = g,
    GAMMA = gamma,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    PVALUE = PVALUE, 
    POWER = POWER
  )
  
  return(result)
  
}

SurvQtestPieceExp_empirical <- function(n1, n2, k1, k2, lamb_group1, lamb_end, tcut, p, cens_type = 'unif', cens_rate = 0.3, alpha = 0.05, late = TRUE) {
  truetimes1_unsort <- rexp(n1, rate = lamb_group1)   # Exponential distribution
  if(cens_type == 'unif'){
    censtimes1_unsort <- runif(n1, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes1_unsort <- rexp(n1, rate = cens_rate)
  }
  
  #####
  status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  data1 <- data.frame("x" = alltimes1_unsort , "event" = status1_unsort) #truetimes0_unsort
  Time1 <- alltimes1_unsort
  ##
  
  #status1_unsort = as.numeric(truetimes1_unsort <= censtimes1_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  #alltimes1_unsort = pmin(truetimes1_unsort, censtimes1_unsort) # Observed times: min between true times and censored times
  
  #Tsort1 = sort(alltimes1_unsort, index.return = TRUE)
  #alltimes1 = Tsort1$x
  #status1 = alltimes1[Tsort1$ix]
  #data1 <- data.frame("x" = alltimes1, "event" = status1) #truetimes0_unsort
  
  data1_PLOT = data.frame(data1)
  data1_PLOT$groups = rep(1, n1)
  
  #quantile_group1 = find_quantile(data1)
  km_fit1 <- survfit(Surv(alltimes1_unsort, status1_unsort) ~ 1)
  quantile_group1 = unname(summary(km_fit1)$table["median"])
  
  lamb_begin = lamb_group1
  
  rate = c(lamb_begin, lamb_end)
  truetimes2_unsort = rsurv(n2, cuts = c(tcut), alpha = rate)
  
  if(cens_type == 'unif'){
    censtimes2_unsort <- runif(n2, min = k2, max = k1+k2)
  }
  else if(cens_type == 'exp'){
    censtimes2_unsort <- rexp(n2, rate = cens_rate)
  }
  
  #status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  #alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  #Tsort2 = sort(alltimes2_unsort, index.return = TRUE)
  #alltimes2 = Tsort2$x
  #status2 = alltimes2[Tsort2$ix]
  
  status2_unsort = as.numeric(truetimes2_unsort <= censtimes2_unsort) # Censoring status: 1 if uncensored, 0 if censored
  
  alltimes2_unsort = pmin(truetimes2_unsort, censtimes2_unsort) # Observed times: min between true times and censored times
  
  data2 <- data.frame("x" = alltimes2_unsort , "event" = status2_unsort) #truetimes0_unsort
  Time2 <- alltimes2_unsort
  ##
  
  #data2 <- data.frame("x" = alltimes2, "event" = status2) #truetimes0_unsort
  
  #Time2 <- alltimes2_unsort
  
  km_fit2 <- survfit(Surv(alltimes2_unsort, status2_unsort) ~ 1)
  quantile_group2 = unname(summary(km_fit2)$table["median"])
  
  #quantile_group2 = find_quantile(data2)
  delta = quantile_group1 - quantile_group2
  
  n1 = length(data1$x)
  n2 = length(data2$x)
  
  n = n1+n2
  mu = n1/n
  
  #fun_integrand = function(rate, x){
  #return(rate*x*exp(rate*x))
  #}
  
  # estimating phi and gamma numerically
  #phi = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group1, rate = lamb_group1)$val # they are the same as computing with the 'true_phi' function
  #gamma = ((1-p)^2)*integrate(fun_integrand, 0, quantile_group2, rate = lamb_group2)$val
  
  # computing explicit phi and gamma
  if(cens_type == 'unif'){
    if(quantile_group1 <k2){
      phi = (1-p)^2*(exp(lamb_group1*quantile_group1) - 1)
      cat('ERROR: Add computation for gamma')
      stop("ERROR")
    }
  }
  
  else if(cens_type == 'exp'){
    phi = (1-p)^2*(lamb_group1/(lamb_group1+cens_rate))*(exp((lamb_group1+cens_rate)*quantile_group1) -1)
    #gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) + (1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin-lamb_end)*tcut))*(exp((lamb_end+cens_rate)*(quantile_group2 - tcut)))
    gamma = (1-p)^2*(lamb_begin/(lamb_begin+cens_rate))*(exp((lamb_begin+cens_rate)*tcut) -1) + (1-p)^2*(lamb_end/(lamb_end+cens_rate))*(exp((lamb_begin+cens_rate)*tcut))*(exp((lamb_end+cens_rate)*(quantile_group2 - tcut)))
    
  }
  
  
  f = lamb_group1*exp(-lamb_group1 * quantile_group1)
  g = g_piecewise(tcut, lamb_begin, lamb_end, quantile_group2)
  
  var_H0 = phi/(mu*(f^2))  + gamma/((1-mu)* g^2)
  
  stat = (delta)/sqrt(var_H0)
  
  stat = sqrt(n)*stat
  
  z <- qnorm(1 - alpha / 2)  # Z-score for the specified alpha level
  STD_ERR <- sqrt(var_H0)
  LOWER <- delta - z * sqrt(var_H0)
  UPPER <- delta + z * sqrt(var_H0)
  PVALUE <- 2*(1-pnorm(abs(stat)) )
  #PVALUE <- 1 - pchisq(stat^2, df = 1) # because the stat is squared, we use chisq instead of gaussian
  POWER <- 1-pnorm(z - sqrt(n/var_H0)*delta ) + pnorm(-z - sqrt(n/var_H0)*delta )
  
  #print(c('stat' = stat, 'std' = sqrt(var_H0), 'phi' = phi, 'f' = f, 'g' = g, 'gamma' = gamma, 'delta' = delta))
  
  result <- data.frame(
    NUM_F = n1,
    NUM_G = n2,
    XI_F = quantile_group1,
    XI_G = quantile_group2,
    delta = delta,
    p = p,
    LAMB_GROUP2 = lamb_end,
    STAT = stat,
    g = g,
    GAMMA = gamma,
    STD_ERR = STD_ERR, 
    LOWER = LOWER, 
    UPPER = UPPER, 
    PVALUE = PVALUE, 
    POWER = POWER
  )
  
  return(result)
  
}

