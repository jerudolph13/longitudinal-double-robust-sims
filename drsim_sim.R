
###################################################################################################
#
# Purpose: Assess the performance of LTMLE under an increasing number of time points and
#           with different sample sizes
#
# Author: Jacqueline Rudolph
#
# Last Update: 01 Oct 2020
#
##################################################################################################

lib <- "~/R/x86_64-pc-linux-gnu-library/4.0"
packages <- c("tidyverse", "ltmle", "SuperLearner", "ranger", "parallel", "gam", "earth")
for (package in packages) {
  library(package, character.only=T, lib.loc=lib)
}


# Pull in command line arguments
args <- commandArgs(trailingOnly=TRUE)

# Define parameters and functions
nsim <- 1000              #Number of simulations
n <- as.numeric(args[1])  #Sample size: {200, 500, 1000, 2000, 5000, 10000}
k <- as.numeric(args[2])  #Number of time points: {1, 5, 10, 20, 50}
yk <- 0.5/k               #Probability of the outcome at any given time point
cores <- detectCores()

'%nin%' = Negate('%in%')


###################################################################################################################################  
#Start simulation

#What is the truth for this scenario?
truth.k <- read.table("drsim_truth.txt", header=TRUE, sep="\t")
true_rd <- truth.k$rd[truth.k$k==k]

res <- data.frame(method=c("glm_default", "glm_user", "SL_default", "SL_user"),
                  sim=rep(NA, 4),
                  r1=rep(NA,4),
                  r0=rep(NA,4),
                  rd=rep(NA,4),
                  se=rep(NA,4))

sim_loop <- function(r) {
  
###################################################################################################################################  
#Generate data
set.seed(123 + r)

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
        XLast.i[j] <- 0
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
                            - 0.4*Z.i[j]  + 0.4*0.5
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
  #Locate first event, after that let Y be 1 (required by LTMLE)
  group_by(ID) %>% 
  mutate(cumy=cumsum(Y), 
         Y=if_else(cumy>=1, 1, 0))

  
###################################################################################################################################  
#Get data in required LTMLE format (wide, with proper variable ordering)

  #Baseline confounders
  baseline <- select(sim.dat[sim.dat$time==1, ], ID, J, K, L)
      base <- paste(names(baseline)[names(baseline) %nin% "ID"], collapse="+")

  #Exposure
  X <- pivot_wider(sim.dat, id_cols="ID", names_from="time", values_from="X", names_prefix="X")
  
  #Time-varying confounder
  Z <- pivot_wider(sim.dat, id_cols="ID", names_from="time", values_from="Z", names_prefix="Z")
  
  #Outcome
  Y <- pivot_wider(sim.dat, id_cols="ID", names_from="time", values_from="Y", names_prefix="Y")
  
  #Now put in correct order
  sim.dat.wide <- baseline
  for (i in 1:k) {
    sim.dat.wide <- merge(sim.dat.wide, Z[ , c(1, i+1)], by="ID")
    sim.dat.wide <- merge(sim.dat.wide, X[ , c(1, i+1)], by="ID")
    sim.dat.wide <- merge(sim.dat.wide, Y[ , c(1, i+1)], by="ID")
  }
  ltmle.dat <- select(sim.dat.wide, -ID)

  
###################################################################################################################################  
#Provide other specifications required for LTMLE

  #Var names
  Anodes <- names(X)[names(X) %nin% "ID"]
  Lnodes <- names(Z)[names(Z) %nin% "ID"]
  Ynodes <- names(Y)[names(Y) %nin% "ID"]
  
  #User-defined formulas
  Aformula <- rep(NA, k)
  for (i in 1:k){
    if (i==1){
      Aformula[i] <- paste(Anodes[i], "~", base, sep="")
    } else {
      Aformula[i] <- paste(Anodes[i], "~", Anodes[i-1], "+", Lnodes[i], "+", Lnodes[i-1], "+", base, sep="")
    }
  }
  
  LYformula <- rep(NA, k*2)
  for (i in 1:k){
    if (i==1){
      #Confounder
      names(LYformula)[2*i-1] <- Lnodes[i]
        LYformula[2*i-1] <- "Q.kplus1 ~ 1"
      #Outcome
      names(LYformula)[2*i] <- Ynodes[i]
        LYformula[2*i] <- paste("Q.kplus1 ~", Anodes[i], "+", base, sep="")
    } else {
      #Confounder
      names(LYformula)[2*i-1] <- Lnodes[i]
        LYformula[2*i-1] <- paste("Q.kplus1 ~", paste(Lnodes[i-1], Anodes[i-1], sep="+"), sep="")
      #Outcome
      names(LYformula)[2*i] <- Ynodes[i]
      LYformula[2*i] <- paste("Q.kplus1 ~", paste(Anodes[i], Anodes[i-1], Lnodes[i], Lnodes[i-1], base, sep="+"), sep="")
    }
  }
  
  #Interventions (always vs. never exposed)
  abar1 <- rep(1, k)
  abar0 <- rep(0, k)
  
  #SuperLearner library
  SL.lib <- c("SL.glm", "SL.gam", "SL.earth")#, "SL.ranger")
  
  
###################################################################################################################################  
#Run LTMLE
  
  #Using GLM models, default specification
  ltmle1 <- ltmle(ltmle.dat, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                  Qform=NULL,
                  gform=NULL,
                  abar=list(treatment=abar1, control=abar0), 
                  survivalOutcome=TRUE, 
                  SL.library=NULL)
  summ <- summary(ltmle1)  
    res$r1[1] <- summ$effect.measures$treatment$estimate
    res$r0[1] <- summ$effect.measures$control$estimate
    res$rd[1] <- summ$effect.measures$ATE$estimate
    res$se[1] <- summ$effect.measures$ATE$std.dev
  
  #Using GLM models, user-defined formula
  ltmle2 <- ltmle(ltmle.dat, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                  Qform=LYformula,
                  gform=Aformula,
                  abar=list(treatment=abar1, control=abar0), 
                  survivalOutcome=TRUE, 
                  SL.library=NULL)
  summ <- summary(ltmle2)  
    res$r1[2] <- summ$effect.measures$treatment$estimate
    res$r0[2] <- summ$effect.measures$control$estimate
    res$rd[2] <- summ$effect.measures$ATE$estimate
    res$se[2] <- summ$effect.measures$ATE$std.dev
  
  #Using SuperLearner, default specification
  ltmle3 <- ltmle(ltmle.dat, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                  Qform=NULL,
                  gform=NULL,
                  abar=list(treatment=abar1, control=abar0), 
                  survivalOutcome=TRUE, 
                  SL.library=SL.lib)
  summ <- summary(ltmle3)  
    res$r1[3] <- summ$effect.measures$treatment$estimate
    res$r0[3] <- summ$effect.measures$control$estimate
    res$rd[3] <- summ$effect.measures$ATE$estimate
    res$se[3] <- summ$effect.measures$ATE$std.dev
    
  #Using SuperLearner, user-defined formulas
  ltmle4 <- ltmle(ltmle.dat, Anodes=Anodes, Lnodes=Lnodes, Ynodes=Ynodes, 
                  Qform=LYformula,
                  gform=Aformula,
                  abar=list(treatment=abar1, control=abar0), 
                  survivalOutcome=TRUE, 
                  SL.library=SL.lib)
  summ <- summary(ltmle4)  
    res$r1[4] <- summ$effect.measures$treatment$estimate
    res$r0[4] <- summ$effect.measures$control$estimate
    res$rd[4] <- summ$effect.measures$ATE$estimate
    res$se[4] <- summ$effect.measures$ATE$std.dev
  
  res$sim <- rep(r, 4)
  return(res)
}   

all.res <- mclapply(1:nsim, function(x) {sim_loop(x)}, mc.cores=cores, mc.set.seed=FALSE)
all.res <- do.call(rbind, all.res)


###################################################################################################################################  
#Summarize and output

#Summary measures
all.res <- all.res %>% 
  mutate(abs.err = abs(rd - true_rd),
         cover = (rd - 1.96*se) <= true_rd & true_rd <= (rd + 1.96*se),
         ciwidth = (rd + 1.96*se) - (rd - 1.96*se))

summ.res <- all.res %>% 
  group_by(method) %>% 
  summarize(avg.rd = mean(rd),
            sd.rd = sd(rd),
            avg.se = mean(se),
            avg.err = mean(abs.err),
            coverage = mean(cover),
            avg.width = mean(ciwidth)) %>% 
  mutate(bias = avg.rd - true_rd,
         mse = bias^2 + sd.rd^2,
         eff = sd.rd/avg.se)

#Output results
filename <- paste("./results/drsim_n-",n,"_k-",k,".txt", sep="")
write.table(all.res, file=filename, sep="\t", row.names=FALSE)
