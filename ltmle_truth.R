
###################################################################################################
#
# Project: Assess the performance of LTMLE under an increasing number of time points and
#           with different sample sizes
#
# Purpose: This program obtains the "true" ATE, using two samples of 1 million (in the first,
#           set exposure to 1; in the second, set to 0)
#
# Author: Jacqueline Rudolph
#
# Last Update: 20 Apr 2020
#
##################################################################################################

setwd("/home/jackie/Documents/drsim/results")

packages <- c("tidyverse", "survival", "parallel")
for (package in packages) {
  update.packages(package, ask=F)
  library(package, character.only=T)
}

'%nin%' = Negate('%in%')


###################################################################################################################################  
#Set up function to get true RD for each of the max time points examined

n <- 500000                 #Sample size
K <- c(1, 5, 10, 20, 50)     #Number of time points

truth <- data.frame(k=NA, rd=NA)

res.k <- function (h) {

  truth$k <- k <- K[h]  #Select number of time points
  yk <- 0.3/k           #Probablity of the outcome at any given time point
  
  
###################################################################################################################################  
#Generate data

create_sim.dat <- function (i, k, exposure) {
    set.seed(n*h + i)
    
    ID <- rep(i, k)
    time <- c(1:k)
    J.i<-K.i<-L.i<-Z.i<-ZLast.i<-X.i<-XLast.i<-Y.i <- rep(0, k)
    
    for (j in 1:k) {
      if (j==1){
        #Baseline confounders
        J.i[j] <- rbinom(1, 1, 0.5)
        K.i[j] <- rbinom(1, 1, 0.5)
        L.i[j] <- rbinom(1, 1, 0.5)
        
        #Time-varying confounder
        ZLast.i[j] <- 0
        Z.i[j] <- rbinom(1, 1, 0.5)
        
        #Set exposure
        XLast.i[j] <- 0
        X.i[j] <- exposure
        
        #Outcome (affected by baseline exposure and confounders)
        p_y <- 1/(1 + exp(-(-log(1/yk - 1) + log(2)*X.i[j]
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        Y.i[j] <- rbinom(1, 1, p_y)
      } else {
        #Baseline confounders
        J.i[j] <- J.i[j-1]
        K.i[j] <- K.i[j-1]
        L.i[j] <- L.i[j-1]
        
        #Time-varying confounder (affected by last Z value and last exposure value)
        ZLast.i[j] <- Z.i[j-1]
        p_z <- 1/(1 + exp(-(-log(1/0.5 - 1) + 0.8*ZLast.i[j]
                            + 0.4*XLast.i[j])))
        Z.i[j] <- rbinom(1, 1, p_z)
        
        #Set exposure
        XLast.i[j] <- X.i[j-1]
        X.i[j] <- exposure
        
        #Outcome (affected by current+last exposures, current+last Z values, baseline confounders)
        p_y <- 1/(1 + exp(-(-log(1/yk - 1) + log(2)*X.i[j] + log(2)*XLast.i[j]
                            -0.4*Z.i[j] - 0.4*ZLast.i[j]
                            + 0.5*J.i[j] - 0.5*0.5
                            + 0.3*K.i[j] - 0.3*0.5
                            - 0.75*L.i[j] + 0.75*0.5)))
        Y.i[j] <- rbinom(1, 1, p_y)
        
      }
    }
  sim.i <- data.frame(ID=ID, time=time, J=J.i, K=K.i, L=L.i, Z=Z.i, ZLast=ZLast.i, X=X.i, XLast=XLast.i, Y=Y.i)
  return(sim.i)
}

#Create exposed sample
cores <- detectCores()
sim.dat1 <- mclapply(1:n, function(x) {create_sim.dat(x, k=k, exposure=1)}, mc.cores=cores, mc.set.seed=FALSE)
sim.dat1 <-do.call(rbind,sim.dat1) %>% 
  #Remove records after the first event (NOTE: this differs from what LTMLE requires)
  group_by(ID) %>% 
  mutate(cumy=cumsum(cumsum(Y))) %>% 
  filter(cumy<=1) %>% 
  select(-cumy)

#Create unexposed sample
sim.dat0 <- mclapply(1:n, function(x) {create_sim.dat(x, k=k, exposure=0)}, mc.cores=cores, mc.set.seed=FALSE)
sim.dat0 <-do.call(rbind,sim.dat0) %>% 
  #Remove records after the first event (NOTE: this differs from what LTMLE requires)
  group_by(ID) %>% 
  mutate(cumy=cumsum(cumsum(Y))) %>% 
  filter(cumy<=1) %>% 
  select(-cumy)

  
###################################################################################################################################  
#Estimate risk in each sample non-parametrically and compare

#First subset to last record
sim.dat.last1 <- sim.dat1 %>% 
  mutate(last=as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)

sim.dat.last0 <- sim.dat0 %>% 
  mutate(last=as.numeric(!duplicated(ID, fromLast=T))) %>% 
  filter(last==1) %>% 
  select(-last)

#Now estimate risk
surv1 <- survfit(Surv(time, Y) ~ 1, data=sim.dat.last1)
  risk1 <- 1 - surv1$surv[length(surv1$surv)]

surv0 <- survfit(Surv(time,Y) ~ 1, data=sim.dat.last0)
  risk0 <- 1 - surv0$surv[length(surv0$surv)]

#Risk difference
truth$rd <- risk1 - risk0

return(truth)

}


###################################################################################################################################  
#Output results

truth.k <- lapply(1:length(K), function(x) {res.k(x)})
  truth.k <- do.call(rbind, truth.k)

write.table(truth.k, file="ltmle_sim_truth.txt", sep="\t", row.names=FALSE)
