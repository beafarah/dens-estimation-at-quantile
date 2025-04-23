rm(list=ls())
require(survival)
require(parallel)
###############################################################################################
# Code for illustrating that there exists a plateau in the MSE around the value of sigma that minimizes the MSE

###############################################################################################
# Generates the MSE for the Exponential distribution
###############################################################################################
# change the sample size in order to see the plateau changing
n <- 1000
lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
censparam <- 0.12
TrueT <- rexp(n,lambda) #true quantile equals 0.5
Cens <- rexp(n,censparam)
status <- TrueT <= Cens #mean(status)# 0.75727#0.75512
ObsT <- pmin(TrueT,Cens)

sigma <- 1
B <- 1e+05 #nb of Gaussian replications
Tb <- rnorm(B,0,sigma)
trueQ <- log(2)/lambda
Yb <- sqrt(n)*(pexp(trueQ+Tb/sqrt(n),lambda)-0.5)
lm_mod <- lm(Yb ~ Tb-1)
summary(lm_mod)

KMfit <- survfit(Surv(ObsT,status) ~ 1)
Qhat <- quantile(KMfit,probs = 0.5)$quantile
KMfun <- stepfun(KMfit$time,c(1,KMfit$surv))
prob <- 1-(KMfun(Qhat+Tb/sqrt(n)))

Yb_hat <- sqrt(n)*(prob-0.5)
lm_mod_hat <- lm(Yb_hat ~ Tb-1)#sort(Tb)
summary(lm_mod_hat)

#setwd(choose.dir())
setwd("C:/Users/beafa/OneDrive/Documents/GitHub/dens-estimation-at-quantile/MSE/n1000")
###################################################################
#setwd("C:/Users/beafa/OneDrive/Documents/Results_Simus_Sig01")
MonteCarlo <- function(n){
  n <- 1000
  lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
  censparam <- 0.12
  TrueT <- rexp(n,lambda) #true quantile equals 0.5
  Cens <- rexp(n,censparam)
  status <- TrueT <= Cens #mean(status)# 0.75727#0.75512
  ObsT <- pmin(TrueT,Cens)
  
  sigma <- 0.1
  B <- 1e+05 #nb of Gaussian replications
  Tb <- rnorm(B,0,sigma)
  trueQ <- log(2)/lambda
  Yb <- sqrt(n)*(pexp(trueQ+Tb/sqrt(n),lambda)-0.5)
  lm_mod <- lm(Yb ~ Tb-1)
  densit_estim <- lm_mod$coefficients
  
  KMfit <- survfit(Surv(ObsT,status)~1)
  Qhat <- quantile(KMfit,probs = 0.5)$quantile
  KMfun <- stepfun(KMfit$time,c(1,KMfit$surv))
  prob <- 1-(KMfun(Qhat+Tb/sqrt(n)))
  #prob <- 1-summary(KMfit,times = (Qhat+Tb/sqrt(n)))$surv
  Yb_hat <- sqrt(n)*(prob-0.5)
  lm_mod_hat <- lm(Yb_hat ~ Tb-1)#sort(Tb)
  densit_estim_hat <- lm_mod_hat$coefficients
  return(list(densit_estim=densit_estim,densit_estim_hat=densit_estim_hat))
}
#MonteCarlo(50)#test
nb_cores <- 8
cl <- makeCluster(nb_cores)
clusterEvalQ(cl, library(survival))
M <- 125 #number of repetitions = M*nb_cores
n <- 200
for (i in 1:M){
  result <- parLapply(cl, X = rep(n, nb_cores), fun = MonteCarlo)
  #result=mclapply(X=rep(n,nb_cores), FUN=MonteCarlo, mc.cores=nb_cores)
  save(result,file=paste("sim",i,".Rdata",sep=""))
}
stopCluster(cl)
#########################################################################

#Extract results
out<-list()
densit_estim <- densit_estim_hat <- rep(NA,M*nb_cores)
for (i in 1:M)
{  
  load(paste("sim",i,".Rdata",sep=""))
  out[[i]]<-result
  res <- unlist(result)
  densit_estim[(8*(i-1)+1):(8*i)] <- res[seq(1,15,by=2)]
  densit_estim_hat[(8*(i-1)+1):(8*i)] <- res[seq(2,16,by=2)]
}

#The true parameter
TrueDens <- lambda/2 
Bias_densit <- mean(densit_estim-TrueDens)
Var_densit <- var(densit_estim)
MSE_densit <- Bias_densit^2+Var_densit

Bias_densit_hat <- mean(densit_estim_hat-TrueDens)
Var_densit_hat <- var(densit_estim_hat)
MSE_densit_hat <- Bias_densit_hat^2+Var_densit_hat

#########################################################################
# Loop on sigma
setwd(choose.dir())
#setwd("/Users/obouaziz/Seafile/Theses/These_Beatriz/Codes_moi/Results_Loop_Sigma")
Sigma <- seq(0.1,15,by=0.1) # sd
M <- 100
densit_estim_hat = rep(NA, M)

#n <- 50
#n <- 200
#n <- 1000

n <- 1000

#B <- 1e+05 #nb of Gaussian replications
B <- 1e+05

lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
censparam <- 0.12

nb_cores = 8
cl <- makeCluster(nb_cores)
clusterEvalQ(cl, library(survival))
clusterExport(cl, varlist = ls())


set.seed(12)
Bias_densit_hat <- Var_densit_hat <- MSE_densit_hat <- rep(NA,(length(Sigma)))
for (k in 1:(length(Sigma))){
  clusterExport(cl, varlist = ls())
  MonteCarlo <- function(n){
    #lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
    #censparam <- 0.48
    TrueT <- rexp(n,lambda) #true quantile equals 0.5
    Cens <- rexp(n,censparam)
    status <- TrueT <= Cens #mean(status)# 0.75727#0.75512
    ObsT <- pmin(TrueT,Cens)
    
    sigma <- Sigma[k]
    
    Tb <- rnorm(B,0,sd = sigma)
    trueQ <- log(2)/lambda
    Yb <- sqrt(n)*(pexp(trueQ+Tb/sqrt(n),lambda)-0.5)
    lm_mod <- lm(Yb ~ Tb-1)
    densit_estim <- lm_mod$coefficients
    
    KMfit <- survfit(Surv(ObsT,status)~1)
    Qhat <- quantile(KMfit,probs = 0.5)$quantile
    KMfun <- stepfun(KMfit$time,c(1,KMfit$surv))
    
    #prob <- 1-(KMfun(trueQ+Tb/sqrt(n))) #####################################################
    prob <- 1-(KMfun(Qhat+Tb/sqrt(n)))
    
    
    #prob <- 1-summary(KMfit,times = (Qhat+Tb/sqrt(n)))$surv
    Yb_hat <- sqrt(n)*(prob-0.5)
    lm_mod_hat <- lm(Yb_hat ~ Tb-1)#sort(Tb)
    densit_estim_hat <- lm_mod_hat$coefficients
    return(list(densit_estim=densit_estim,densit_estim_hat=densit_estim_hat))
  }
  #for (i in 1:M){
    #resulti=lapply(X=rep(n,M), FUN=MonteCarlo)
    #result=parLapply(cl, X = rep(n, nb_cores), fun = MonteCarlo)
    #save(result,file=paste("sim",i,".Rdata",sep=""))
  #}
  set.seed(k)
  resulti=parLapply(cl, X = rep(n, M), fun = MonteCarlo)
  #resulti=lapply(X=rep(n,M), FUN=MonteCarlo)
  
  #out<-list()
  #densit_estim_hat <- rep(NA,M*nb_cores)
  for (i in 1:M)
  {  
    densit_estim_hat[i] = resulti[[i]]$densit_estim_hat
    #load(paste("sim",i,".Rdata",sep=""))
    #out[[i]]<-result
    #res <- unlist(result)
    #densit_estim_hat[(8*(i-1)+1):(8*i)] <- res[seq(2,16,by=2)]
  }
  TrueDens <- lambda/2
  Bias_densit_hat[k] <- mean(densit_estim_hat)-TrueDens
  Var_densit_hat[k] <- var(densit_estim_hat)
  MSE_densit_hat[k] <- (Bias_densit_hat[k])^2+Var_densit_hat[k]
  print(k)
}
save(Bias_densit_hat,Var_densit_hat,MSE_densit_hat,file="results_loop_sigma.RData")
stopCluster(cl)

#file_path <- file.choose()
#load("results_loop_sigma.RData")

plot(Sigma,Bias_densit_hat,type="p")
plot(Sigma,Var_densit_hat,type="p")
plot(Sigma,MSE_densit_hat,type="p")

plot(Sigma,MSE_densit_hat,type="b", main = paste("n = ", n ), ylab = 'MSE', xlab = expression(sigma) ) 

# sigma that minimizes MSE
index_min = which.min(MSE_densit_hat)
cat(paste("Sigma that minimizes MSE is: "), Sigma[index_min])

#######
#load("C:/Users/beafa/OneDrive/Documents/GitHub/these/density estimation/Results_Loop_Sigma/results_loop_sigma.RData")
Sigma <- seq(0.1,15,by=0.1)
plot(Sigma, MSE_densit_hat, type="b", 
     main = paste("n =", n), 
     xlab = expression(sigma), 
     ylab = "MSE")

load("C:/Users/beafa/OneDrive/Documents/GitHub/these/density estimation/Results_Loop_Sigma_n1000/results_loop_sigma.RData")
plot(Sigma,MSE_densit_hat, type="p", main='n=1000')




##########################################################################################
sigmas_lst = seq(0.1, 15, by = 0.1)

n <- 50
lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
censparam <- 0.48
B <- 1e+05

first_lst = c()
second_lst = c()
third_lst = c()

total_lst = c()
for(sigma in sigmas_lst){
  TrueT <- rexp(n,lambda) #true quantile equals 0.5
  Cens <- rexp(n,censparam)
  status <- TrueT <= Cens #mean(status)# 0.75727#0.75512
  ObsT <- pmin(TrueT,Cens)
  
  KMfit <- survfit(Surv(ObsT,status) ~ 1)
  Qhat <- quantile(KMfit,probs = 0.5)$quantile # estimation of the quantile

  Z = rnorm(1, 0, 1)
  first_term = exp(-lambda*Qhat)*sqrt((lambda/(lambda+censparam)))*exp((lambda+censparam)*Qhat/2)*sqrt(((lambda+censparam)/sqrt(n)))*Z*( 2*gamma(5/4)/(sqrt(sigma)*sqrt(2*pi)*2^(-1/4))  + (lambda+censparam)*(2*sqrt(sigma) * gamma(7/4))/(4*sqrt(n)*sqrt(2*pi)*2^(-3/4)))
  
  second_term = -( (lambda/(lambda+censparam))^(1/2)*(exp(Qhat*(lambda+censparam)) -1 )^(1/2)*Z*(  (lambda+censparam)/sqrt(n) * (exp(Qhat*(lambda+censparam)))/(2*(exp(Qhat*(lambda+censparam)) - 1) ) *( ( (3*sigma^2)*(-lambda^2*exp(-lambda*Qhat)) )/(2*n)  ) + ((lambda*exp(-lambda*Qhat))/(sigma^2*sqrt(n) ) + 3*lambda^3*exp(-lambda*Qhat)/(6*n^(3/2)) ) ) )  
  third_term = lambda*exp(-lambda*Qhat) + (3*sigma^2*lambda^3*exp(-lambda*Qhat))/(6*n^(3/2)) 
  
  first_lst = c(first_lst, first_term)
  second_lst = c(second_lst, second_term)
  third_lst = c(third_lst, third_term)

  total_lst = c(total_lst, first_term + second_term + third_term)
}


plot(sigmas_lst, total_lst)
index_min = which.min(total_lst)

sigmas_lst[index_min]

plot(sigmas_lst, second_lst)


##########################################
# Multiple repetitions of the code

sigmas_lst = seq(0.1, 15, by = 0.1)

n <- 1000
lambda <- 1.5 # median = log(2)/lambda = 0.4620981 #dexp(log(2)/lambda,1.5) =0.75 #2*log(2)
censparam <- 0.48
B <- 1e+05

M = 500
it = 0

mean_first = numeric(length(sigmas_lst))
mean_second = numeric(length(sigmas_lst))
mean_third = numeric(length(sigmas_lst))
mean_total = numeric(length(sigmas_lst))
while(it < M){
  print(it)
  it = it+1
  
  for(i in seq_along(sigmas_lst)){
    sigma = sigmas_lst[i]
    TrueT <- rexp(n,lambda) #true quantile equals 0.5
    Cens <- rexp(n,censparam)
    status <- TrueT <= Cens #mean(status)# 0.75727#0.75512
    ObsT <- pmin(TrueT,Cens)
    
    KMfit <- survfit(Surv(ObsT,status) ~ 1)
    Qhat <- quantile(KMfit,probs = 0.5)$quantile # estimation of the quantile
    
    
    ############# WITH TRUE QUANTILE
    trueQ <- log(2)/lambda
    #Qhat = trueQ
    
    Z = rnorm(1, 0, 1)
    
    first_term = exp(-lambda*Qhat)*sqrt((lambda/(lambda+censparam)))*exp((lambda+censparam)*Qhat/2)*sqrt(((lambda+censparam)/sqrt(n)))*Z*( 2*gamma(5/4)/(sqrt(sigma)*sqrt(2*pi)*2^(-1/4))  + (lambda+censparam)*(2*sqrt(sigma) * gamma(7/4))/(4*sqrt(n)*sqrt(2*pi)*2^(-3/4)))
    
    second_term = -( (lambda/(lambda+censparam))^(1/2)*(exp(Qhat*(lambda+censparam)) -1 )^(1/2)*Z*(  (lambda+censparam)/sqrt(n) * (exp(Qhat*(lambda+censparam)))/(2*(exp(Qhat*(lambda+censparam)) - 1) ) *( ( (3*sigma^2)*(-lambda^2*exp(-lambda*Qhat)) )/(2*n)  ) + ((lambda*exp(-lambda*Qhat))/(sigma^2*sqrt(n) ) + 3*lambda^3*exp(-lambda*Qhat)/(6*n^(3/2)) ) ) )  
    third_term = lambda*exp(-lambda*Qhat) + (3*sigma^2*lambda^3*exp(-lambda*Qhat))/(6*n^(3/2)) 
    
    mean_first[i] = mean_first[i] + first_term
    mean_second[i] = mean_second[i] + second_term
    mean_third[i] = mean_third[i] + third_term

    
    mean_total[i] = mean_total[i] + first_term + second_term +third_term
    
  }

}

mean_first = mean_first/M
mean_second = mean_second/M
mean_third = mean_third/M
mean_total = mean_total/M

#plot(sigmas_lst, mean_first)
plot(sigmas_lst, mean_total, main = 'n = 1000')

index_min = which.min(mean_total)

sigmas_lst[index_min]
