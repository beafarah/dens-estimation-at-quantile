library(survival)
library(survminer)
require(survival)
library(dplyr)
library(survminer)
library(latex2exp)
library(viridisLite)
library(ggsci) 
library(RColorBrewer)
#rm(list = ls())
setwd("C:/Users/beafa/OneDrive/Documents/GitHub/dens-estimation-at-quantile/src")
source("code_source.R")
source("rsurv.R")
##########################################################################################
# Planning a clinical trial

##########################################################################################
#1) Exponential survival and exponential censoring
alpha = 0.05
lamb_group1 = 1.5
cens_rate = 0.48
#p = 0.1
p = 0.5
n1 = 100
n2 = 100
n = n1 + n2
mu = n1/n

delta = 0.2
k1 = 0 # we dont use these, they are needed only for uniform censoring
k2 = 0
quantile_group1 = -log(1-p)/lamb_group1
quantile_group2 = -delta + quantile_group1
lamb_group2 = -log(1-p)/quantile_group2

# true curves
ts = seq(0, 5, by = 0.01)
#plot(ts, exp(-lamb_group1*ts), type = 'l', lwd = 2, col = 'red', ylab = 'Survival', xlab = 'Time')
plot(ts, exp(-lamb_group1*ts), 
     xlim = c(0,3),
     type = 'l', lwd = 2, col = 'red', ylab = 'Survival', xlab = 'Time',
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5,
     main = TeX(sprintf(r'($\Delta$ = %.2f)', delta)))
lines( ts, exp(-lamb_group2*ts), lwd = 2, lty=2, col = 'blue')
abline(h = 1-p, lty = 2)

######
# delta = 0.2 at p = 0.5
ns = seq(50, 1000, by = 10)
pwrs1 = c()
for(ni in ns){
  res_power1 = SurvQtestExp_new(ni, ni, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
  pwrs1 = c(pwrs1, res_power1)
}

plot(ns, pwrs1, ylim = c(0, 1), type = 'l', lwd = 2, 
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5,
     col = 'darkorange', ylab = "Power", xlab = "Sample size", main = TeX(sprintf(r'($\Delta$ = %.2f)', delta)))
abline(h = 1, lty = 2)


######
# power for different sample sizes: verifying that the power for the explicit test behaves as expected for many deltas 
n_test1 = 50
n_test2 = 100
n_test3 = 300
n_test4 = 500
n_test5 = 1000
#deltas = seq(0, 0.25, by = 0.01)
deltas = seq(0, 1, by = 0.01)
power_deltas1 = c()
power_deltas2 = c()
power_deltas3 = c()
power_deltas4 = c()
power_deltas5 = c()
for(delta in deltas){
  print(delta)
  
  res_explicit1 = SurvQtestExp_new(n_test1, n_test1, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  res_explicit2 = SurvQtestExp_new(n_test2, n_test2, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate) 
  res_explicit3 = SurvQtestExp_new(n_test3, n_test3, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate) 
  res_explicit4 = SurvQtestExp_new(n_test4, n_test4, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate) 
  res_explicit5 = SurvQtestExp_new(n_test5, n_test5, k1, k2, lamb_group1, p, delta, alpha, cens_type = 'exp', cens_rate = cens_rate) 
  
  
  res_power1 = res_explicit1$POWER
  res_power2 = res_explicit2$POWER
  res_power3 = res_explicit3$POWER
  res_power4 = res_explicit4$POWER
  res_power5 = res_explicit5$POWER
  
  power_deltas1 = c(power_deltas1, res_power1)
  power_deltas2 = c(power_deltas2, res_power2)
  power_deltas3 = c(power_deltas3, res_power3)
  power_deltas4 = c(power_deltas4, res_power4)
  power_deltas5 = c(power_deltas5, res_power5)
}

colors <- viridis(5, option = "E") 
colors[5] <- rgb(255, 220, 100, maxColorValue = 255) 
#colors <- pal_jco("default")(5)
#colors <- brewer.pal(5, "Paired")
#colors <- gray(seq(0, 0.8, length.out = 5))
line_types <- c(2, 3, 4, 5, 6)
line_widths = c(4, 4, 3, 2, 2)

plot(deltas, power_deltas1,
     cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1.2,
     type = 'l', lwd = line_widths[1], col = colors[5], 
     lty = line_types[1], 
     main = "", xlab = TeX(r'($\Delta$)'), ylab = "Power", ylim = c(0, 1.0))
lines(deltas,power_deltas2, lwd = line_widths[2], col = colors[4], lty = line_types[2])
lines(deltas,power_deltas3, lwd = line_widths[3], col = colors[3], lty = line_types[3])
lines(deltas,power_deltas4, lwd = line_widths[4], col = colors[2], lty = line_types[4])
lines(deltas,power_deltas5, lwd = line_widths[5], col = colors[1], lty = line_types[5])

legend("bottomright", 
       legend = c(n_test1, n_test2, n_test3, n_test4,n_test5), 
       title=TeX(r'($n_i$)'),
       col = rev(colors), lwd = line_widths, 
       lty = line_types,
       text.font = 1.2,  # Bold text
       cex = 1.2,  # Increase text size
       bty = "o",  # Box around legend
       inset = c(0.01, 0.01))  # Move legend slightly inside


#######
# analytical power
n1 = 100
n2 = 100
delta0 = 0
delta1 = 0.1
delta2 = 0.2
p = 0.5

#analytical
res_power0 = SurvQtestExp_new(n1, n1, k1, k2, lamb_group1, p, delta0, alpha, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
res_power1 = SurvQtestExp_new(n1, n1, k1, k2, lamb_group1, p, delta1, alpha, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
res_power2 = SurvQtestExp_new(n1, n1, k1, k2, lamb_group1, p, delta2, alpha, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time


##########################################################################################
#empirical:

n_its = 10000
#delta = 0
#lamb_group2 = 1.5

# under H0: 
lamb_group2 = 1.5

n1 = 100
type1_err = 0
it = 0
all_simus = data.frame()
pval_lst = c()
stat_lst = c()
while(it < n_its){
  res = SurvQtestExp_empirical(n1, n1, k1, k2, lamb_group1, lamb_group2, p, alpha, cens_type = 'exp', cens_rate = cens_rate)
  all_simus <-rbind(all_simus, res)
  pval = res$PVALUE
  pval_lst = c(pval_lst, pval)
  stat_lst = c(stat_lst, res$STAT)
  if(pval < alpha){
    type1_err = type1_err + 1
  }
  it = it+1
}

pval_unl <- unlist(all_simus$PVALUE, use.names = FALSE)
#hist(pval_unl, main =" ", xlab = "pvalue")

numerator = sum(pval_unl<alpha)
denominator = length(pval_unl)
rejection_rate = numerator/denominator
cat("Percentage of times that wrongfully rejected H0: ", rejection_rate)

######################################################################
# under delta = 0.1 and delta = 0.2
n_its = 10000
delta_pwr = 0.2
quantile_group1 = -log(1-p)/lamb_group1
lamb_group2 = -log(1-p)/(quantile_group1 - delta_pwr)

n1 = 100
type1_err = 0
it = 0
all_simus = data.frame()
pval_lst = c()
stat_lst = c()
while(it < n_its){
  res = SurvQtestExp_empirical(n1, n1, k1, k2, lamb_group1, lamb_group2, p, alpha, cens_type = 'exp', cens_rate = cens_rate)
  all_simus <-rbind(all_simus, res)
  pval = res$PVALUE
  pval_lst = c(pval_lst, pval)
  stat_lst = c(stat_lst, res$STAT)
  if(pval < alpha){
    type1_err = type1_err + 1
  }
  it = it+1
}

pval_unl <- unlist(all_simus$PVALUE, use.names = FALSE)
#hist(pval_unl, main =" ", xlab = "pvalue")

numerator = sum(pval_unl<alpha)
denominator = length(pval_unl)
rejection_rate = numerator/denominator
cat("Empirical power: ", rejection_rate)


##########################################################################################
#1) Piecewise survival and exponential censoring
alpha = 0.05
lamb_group1 = 1.5

delta = 0.1
k1 = 0 # we dont use these, they are needed only for uniform censoring
k2 = 0
tcut = 0.2
p = 0.5
cens_rate = 0.48
#p = 0.1
p = 0.5
n1 = 100
n2 = 100
n = n1 + n2
mu = n1/n
late = TRUE

dat = data_proportional_piece_Exp(n1, n1, k1, k2, lamb_group1, tcut, p, delta)
lamb_end = dat$lamb_group2

if(dat$quantile_group1 - tcut < delta){
  print("Specify another value for delta or tcut!")
}

# true curves
ts = seq(0, 5, by = 0.01)
plot(ts, exp(-lamb_group1*ts), type = 'l', lwd = 2, 
     xlim = c(0,3),
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5,
     col = 'red', ylab = 'Survival', xlab = 'Time', main = TeX(sprintf(r'($\Delta$ = %.2f)', delta)))
lines( ts, sapply(ts, function(x) S_piecewise(tcut, lamb_group1, lamb_end, x)), lwd = 2, lty=2, col = 'blue')
abline(h = 1-p, lty = 2)
#legend(x = "topright",legend = c('Group 1', 'Group 2'), 
#       col = c("red", "blue"), lty = c(1,2) )


# we need quantile_group1 - tcut > delta, or quantile_group2 > tcut
##############################################################################

# delta = 0.2 at p = 0.5
ns = seq(50, 1000, by = 10)
pwrs1_piece = c()
for(ni in ns){
  res_power1_piece = SurvQtestPieceExp_new(ni, ni, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
  pwrs1_piece = c(pwrs1_piece, res_power1_piece)
}

plot(ns, pwrs1_piece, ylim = c(0, 1), 
     cex.lab=1.5, cex.axis=1, cex.main=1.5, cex.sub=1.5,
     type = 'l', lwd = 2, col = 'darkorange', ylab = "Power", xlab = "Sample size", main = TeX(sprintf(r'($\Delta$ = %.2f)', delta)))
abline(h = 1, lty = 2)

######

# power for different sample sizes: verifying that the power for the explicit test behaves as expected for many deltas 
n_test1 = 50
n_test2 = 100
n_test3 = 300
n_test4 = 500
n_test5 = 1000
deltas = seq(0, 0.4, by = 0.01)
power_deltas1 = c()
power_deltas2 = c()
power_deltas3 = c()
power_deltas4 = c()
power_deltas5 = c()
last_delta = tail(deltas, n = 1)
for(delta in deltas){
  print(delta)
  dat = data_proportional_piece_Exp(n_test1, n_test1, k1, k2, lamb_group1, tcut, p, delta)
  
  if(dat$quantile_group1 - tcut < delta){
    last_delta = delta
    break
  }
  
  res_explicit1 = SurvQtestPieceExp_new(n_test1, n_test1, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  res_explicit2 = SurvQtestPieceExp_new(n_test2, n_test2, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  res_explicit3 = SurvQtestPieceExp_new(n_test3, n_test3, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  res_explicit4 = SurvQtestPieceExp_new(n_test4, n_test4, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  res_explicit5 = SurvQtestPieceExp_new(n_test5, n_test5, k1, k2, lamb_group1, tcut, p, delta, alpha, late, cens_type = 'exp', cens_rate = cens_rate) # generates data each time
  
  
  res_power1 = res_explicit1$POWER
  res_power2 = res_explicit2$POWER
  res_power3 = res_explicit3$POWER
  res_power4 = res_explicit4$POWER
  res_power5 = res_explicit5$POWER
  
  power_deltas1 = c(power_deltas1, res_power1)
  power_deltas2 = c(power_deltas2, res_power2)
  power_deltas3 = c(power_deltas3, res_power3)
  power_deltas4 = c(power_deltas4, res_power4)
  power_deltas5 = c(power_deltas5, res_power5)
}


colors <- viridis(5, option = "E") 
colors[5] <- rgb(255, 220, 100, maxColorValue = 255) 
#colors <- pal_jco("default")(5)
#colors <- brewer.pal(5, "Paired")
#colors <- gray(seq(0, 0.8, length.out = 5))
line_types <- c(2, 3, 4, 5, 6)
line_widths = c(4, 4, 3, 2, 2)

plot(deltas[deltas < last_delta], power_deltas1,
     cex.lab=1.5, cex.axis=1, cex.main=1, cex.sub=1.2,
     type = 'l', lwd = line_widths[1], col = colors[5], 
     lty = line_types[1], 
     main = "", xlab = TeX(r'($\Delta$)'), ylab = "Power", ylim = c(0, 1.0))
lines(deltas[deltas < last_delta],power_deltas2, lwd = line_widths[2], col = colors[4], lty = line_types[2])
lines(deltas[deltas < last_delta],power_deltas3, lwd = line_widths[3], col = colors[3], lty = line_types[3])
lines(deltas[deltas < last_delta],power_deltas4, lwd = line_widths[4], col = colors[2], lty = line_types[4])
lines(deltas[deltas < last_delta],power_deltas5, lwd = line_widths[5], col = colors[1], lty = line_types[5])

legend("bottomright", 
       legend = c(n_test1, n_test2, n_test3, n_test4,n_test5), 
       title=TeX(r'($n_i$)'),
       col = rev(colors), lwd = line_widths, 
       lty = line_types,
       text.font = 1.2,  # Bold text
       cex = 1.2,  # Increase text size
       bty = "o",  # Box around legend
       inset = c(0.01, 0.01))  # Move legend slightly inside
#######
# analytical power
n1 = 100
n2 = 100
delta0 = 0
delta1 = 0.1
delta2 = 0.2
p = 0.5

#analytical
res_power0_piece = SurvQtestPieceExp_new(n1, n1, k1, k2, lamb_group1, tcut, p, delta0, alpha, late, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
res_power1_piece = SurvQtestPieceExp_new(n1, n1, k1, k2, lamb_group1, tcut, p, delta1, alpha, late, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time
res_power2_piece = SurvQtestPieceExp_new(n1, n1, k1, k2, lamb_group1, tcut, p, delta2, alpha, late, cens_type = 'exp', cens_rate = cens_rate)$POWER # generates data each time


#######################
#empirical:

n_its = 10000
#delta = 0

###
# under H0: 
lamb_group2 = 1.5

n1 = 100
type1_err = 0
it = 0
all_simus = data.frame()
pval_lst = c()
stat_lst = c()
while(it < n_its){
  res = SurvQtestPieceExp_empirical(n1, n1, k1, k2, lamb_group1, lamb_group2, tcut, p, cens_type = 'exp', cens_rate = cens_rate, alpha = alpha, late = TRUE)
  all_simus <-rbind(all_simus, res)
  pval = res$PVALUE
  pval_lst = c(pval_lst, pval)
  stat_lst = c(stat_lst, res$STAT)
  if(pval < alpha){
    type1_err = type1_err + 1
  }
  it = it+1
}

pval_unl <- unlist(all_simus$PVALUE, use.names = FALSE)
hist(pval_unl, main =" ", xlab = "pvalue")

numerator = sum(pval_unl<alpha)
denominator = length(pval_unl)
rejection_rate = numerator/denominator
cat("Percentage of times that wrongfully rejected H0: ", rejection_rate)


### for delta = 0.1 and delta = 0.2
delta_pwr = 0.2
dat = data_proportional_piece_Exp(n1, n1, k1, k2, lamb_group1, tcut, p, delta_pwr)
lamb_end = dat$lamb_group2


n1 = 100
type1_err = 0
it = 0
all_simus = data.frame()
pval_lst = c()
stat_lst = c()
while(it < n_its){
  res = SurvQtestPieceExp_empirical(n1, n1, k1, k2, lamb_group1, lamb_end, tcut, p, cens_type = 'exp', cens_rate = cens_rate, alpha = alpha, late = TRUE)
  all_simus <-rbind(all_simus, res)
  pval = res$PVALUE
  pval_lst = c(pval_lst, pval)
  stat_lst = c(stat_lst, res$STAT)
  if(pval < alpha){
    type1_err = type1_err + 1
  }
  it = it+1
}

pval_unl <- unlist(all_simus$PVALUE, use.names = FALSE)
#hist(pval_unl, main =" ", xlab = "pvalue")

numerator = sum(pval_unl<alpha)
denominator = length(pval_unl)
rejection_rate = numerator/denominator
cat("Empirical power: ", rejection_rate)
