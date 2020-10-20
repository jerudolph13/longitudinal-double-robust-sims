
###################################################################################################
#
# Purpose: Re-run analysis with IPW to compare against LTMLE
#
# Author: Jacqueline Rudolph
#
# Last Update: 20 Oct 2020
#
##################################################################################################

packages <- c("tidyverse", "survival", "parallel")
for (package in packages) {
  library(package, character.only=T))
}

# Define parameters and functions
nsim <- 1000              #Number of simulations
nboot <- 200              #Number of bootstrap resamples
n <- 1000                 #Sample size: {200, 500, 1000, 2000, 5000, 10000}
k <- 5                    #Number of time points: {1, 5, 10, 20, 50}
yk <- 0.5/k               #Probability of the outcome at any given time point
cores <- detectCores()

'%nin%' = Negate('%in%')


###################################################################################################################################  
#Start simulation

#What is the truth for this scenario?
truth.k <- read.table("../results/drsim_truth.txt", header=TRUE, sep="\t")
true_rd <- truth.k$rd[truth.k$k==k]

sim.res <- data.frame(sim=NA, rd=NA, risk1=NA, risk0=NA, se=NA, avg.wt=NA, min.wt=NA, max.wt=NA)

sim_loop <- function(r) {
  
###################################################################################################################################  
#Generate data
set.seed(123 + r)
sim.res$sim <- r

create_sim.dat <- function (i, k) {
    ID <- rep(i, k)
    time <- c(1:k)
    J.i<-K.i<-L.i<-Z.i<-ZLast.i<-X.i<-XLast.i<-Y.i <- rep(0, k)
    
    for (j in 1:k) {
      if (j==1){
        #Baseline confounders
        J.i[j] <- rbinom(1, 1, 0.5)
        K.i[j] <- rbinom(1, 1, 0.5)
        L.i[j] <- rbinom(1, 1, 0.5)
        
        #Lagged variables
        ZLast.i[j] <- 0
        XLast.i[j] <- 0
        
        #Time-varying confounder
        Z.i[j] <- rbinom(1, 1, 0.5)
        
        #Exposure (affected by baseline confounders)
        p_x <- 1/(1 + exp(-(-log(1/0.5 - 1)
                            - 0.4*Z.i[j] + 0.4*0.5
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        X.i[j] <- rbinom(1, 1, p_x)
        
        #Outcome (affected by baseline exposure and confounders)
        p_y <- 1/(1 + exp(-(-log(1/yk - 1) 
                            + log(2)*X.i[j]
                            - 0.4*Z.i[j] + 0.4*0.5
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        Y.i[j] <- rbinom(1, 1, p_y)
      } else {
        #Baseline confounders
        J.i[j] <- J.i[j-1]
        K.i[j] <- K.i[j-1]
        L.i[j] <- L.i[j-1]
        
        #Lagged variables
        ZLast.i[j] <- Z.i[j-1]
        XLast.i[j] <- X.i[j-1]
        
        #Time-varying confounder (affected by last Z value and last exposure value)
        p_z <- 1/(1 + exp(-(-log(1/0.5 - 1) 
                            + 0.8*ZLast.i[j]
                            + 0.4*XLast.i[j])))
        Z.i[j] <- rbinom(1, 1, p_z)
        
        #Exposure (affected by last exposure value, current+last Z value, baseline confounders)
        p_x <- 1/(1 + exp(-(-log(1/0.5 - 1) 
                            + 0.8*XLast.i[j]
                            - 0.4*Z.i[j] + 0.4*0.5
                            - 0.4*ZLast.i[j]
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        X.i[j] <- rbinom(1, 1, p_x)
        
        #Outcome (affected by current+last exposures, current+last Z values, baseline confounders)
        p_y <- 1/(1 + exp(-(-log(1/yk - 1) 
                            + log(2)*X.i[j] 
                            + log(2)*XLast.i[j]
                            - 0.4*Z.i[j] + 0.4*0.5
                            - 0.4*ZLast.i[j]
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        Y.i[j] <- rbinom(1, 1, p_y)
        
      }
    }
  sim.i <- data.frame(ID=ID, time=time, J=J.i, K=K.i, L=L.i, Z=Z.i, ZLast=ZLast.i, X=X.i, XLast=XLast.i, Y=Y.i)
  return(sim.i)
}

sim.dat <- lapply(1:n, function(x) {create_sim.dat(x, k=k)})
sim.dat <-do.call(rbind,sim.dat) %>% 
  #Remove records after the first event (NOTE: this differs from what LTMLE requires)
  group_by(ID) %>% 
  mutate(cumy=cumsum(cumsum(Y))) %>% 
  filter(cumy<=1) %>% 
  select(-cumy)


###################################################################################################################################  
#Bootstrap sample

boot.res <- data.frame(rd=NA, risk1=NA, risk0=NA, boot=NA, avg.wt=NA, min.wt=NA, max.wt=NA)

bootrep <- function(b) {
  boot.res$boot <- b
  
  #Sample with replacement from the individuals
  firstobs <- sim.dat[sim.dat$time == 1, ]
  samp <- table(firstobs[sample(1:nrow(firstobs),nrow(firstobs),replace=T), (names(sim.dat) == "ID")])
  
  #Grab other records for selected individuals
  boot <- NULL
  if(b==0){
    boot<-sim.dat %>% rename(bid=ID)
  } else{
    for(zzz in 1:max(samp)){ 
      cc <- sim.dat[sim.dat$ID %in% names(samp[samp %in% c(zzz:max(samp))]),]
      cc$bid<-paste0(cc$ID,zzz)
      boot <- rbind(boot, cc)
    }}
  
  
###################################################################################################################################  
#Estimate IP weights (within bootstrap resample)
  
  if (k==1) {
    #Denominator of weights
    ps <- glm(X ~ Z + J + K + L, family=binomial(link="logit"), data=boot)$fitted.values
    boot$denominator <- boot$X*ps + (1-boot$X)*(1-ps)

    #Numerator of weights
    ps <- glm(X ~ 1, family=binomial(link="logit"), data=boot)$fitted.values
    boot$numerator <- boot$X*ps + (1-boot$X)*(1-ps)
  } else {
    #Denominator of weights
    ps <- glm(X ~ XLast:(time>1) + Z + ZLast:(time>1) + J + K + L + as.factor(time), 
              family=binomial(link="logit"), data=boot)$fitted.values
    boot$denominator <- boot$X*ps + (1-boot$X)*(1-ps)
    
    #Numerator of weights
    ps <- glm(X ~ XLast:(time>1) + as.factor(time), family=binomial(link="logit"), data=boot)$fitted.values
    boot$numerator <- boot$X*ps + (1-boot$X)*(1-ps)
  }

  #Divide numerator/denominator and  multiply across time to get weights
  boot <- boot %>% 
    group_by(bid) %>%  
    mutate(wt=cumprod(numerator/denominator)) %>% 
    ungroup(bid) 
    boot.res$avg.wt <- mean(boot$wt)
    boot.res$min.wt <- min(boot$wt)
    boot.res$max.wt <- max(boot$wt)
  
  
###################################################################################################################################  
#Estimate weighted RD at the end of follow-up (within bootstrap)
    #This code was adapted from Hernan and Robins
  
  #Model hazard of the outcome, weighted by IP weights
  if (k==1) {
    haz.mod <- glm(Y ~ X, family=binomial(link="logit"), weight=wt, data=boot)
  } else {
    haz.mod <- glm(Y ~ X + XLast:(time>1) + as.factor(time), family=binomial(link="logit"), weight=wt, data=boot)
  }

  ipw0 <- data.frame(cbind(seq(1, k), 0, 0))
    colnames(ipw0) <- c("time", "X", "XLast")
    
  ipw1 <- data.frame(cbind(seq(1, k), 1, 1))
    colnames(ipw1) <- c("time", "X", "XLast")

  #Assignment of estimated hazard to each time point
  ipw0$p.event0 <- predict(haz.mod, ipw0, type="response")
  ipw1$p.event1 <- predict(haz.mod, ipw1, type="response")

  #Computation of risk for each time point
  ipw0$risk0 <- 1 - cumprod(1 - ipw0$p.event0)
  ipw1$risk1 <- 1 - cumprod(1 - ipw1$p.event1)
  boot.res$rd <- max(ipw1$risk1) - max(ipw0$risk0)
  boot.res$risk1 <- max(ipw1$risk1)
  boot.res$risk0 <- max(ipw0$risk0)
  
  return(boot.res)
}


###################################################################################################################################  
#Summarize across bootstraps

  all.boot <- lapply(0:nboot, function(tt) {bootrep(tt)})
  all.boot <- do.call(rbind, all.boot)

  #Pull out results from original sample to get point estimate
  boot0 <- filter(all.boot, boot == 0)
  sim.res$rd <- boot0$rd
  sim.res$risk1 <- boot0$risk1
  sim.res$risk0 <- boot0$risk0

  #Estimate standard error
  boot.summ <- all.boot %>% 
    summarize(b.sd.rd = sd(rd), b.avg.wt=mean(avg.wt), b.min.wt=min(min.wt), b.max.wt=max(max.wt))
  sim.res$se <- boot.summ$b.sd.rd
  sim.res$avg.wt <- boot.summ$b.avg.wt
  sim.res$min.wt <- boot.summ$b.min.wt
  sim.res$max.wt <- boot.summ$b.max.wt
  
  return(sim.res)
}   

all.res <- mclapply(1:nsim, function(x) {sim_loop(x)}, mc.cores=cores, mc.set.seed=FALSE)
all.res <- do.call(rbind, all.res)


###################################################################################################################################  
#Summarize and output

#Summary measures
all.res <- all.res %>% 
  mutate(cover = (rd - 1.96*se) <= true_rd & true_rd <= (rd + 1.96*se),
         ciwidth = (rd + 1.96*se) - (rd - 1.96*se))

summ.res <- all.res %>% 
  summarize(avg.rd = mean(rd),
            avg.r1 = mean(risk1),
            avg.r0 = mean(risk0),
            sd.rd = sd(rd),
            avg.se = mean(se),
            coverage = mean(cover),
            avg.width = mean(ciwidth)) %>% 
  mutate(bias = avg.rd - true_rd,
         mse = bias^2 + sd.rd^2,
         eff = sd.rd/avg.se)

#Output results
filename <- paste("../results/drsim_ipw_n-",n,"_k-", k,".txt", sep="")
write.table(all.res, file=filename, sep="\t", row.names=FALSE)
 
