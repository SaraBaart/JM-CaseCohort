####################
## Load libraries ##
####################
library(JMbayes)
library(jagsUI)
library(CaseCohort)
library(splines)
library(MASS)
library(reshape2)
library(dplyr)
library(latticeExtra)
library(beepr)
############################
## Obtain Formulas needed ##
############################

# Set working directory
setwd("V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03")
# Simulate data
source("Functions/Simulate-15082017.R")
# Prepare Data for output
source("Functions/PrepJAGS-new.R")
# Fit model with JAGS
#source("Functions/Model Splines.R")
# Extract Extra elements for prediction
source("Functions/ExtractExtra.R")
# Write JAGS as JM object for predictions
source("Functions/JAGSTOJM.R")
# Params
params <- c("betas", "alphas", "gammas", "inv.D", "b", "Bs.gammas", "tau")
# Get CaseCohort data
source("Functions/CaseCohort.R")


##############################
## Test simulation function ##
##############################

test <- Simulate(n = 2000, mean.Cens = 2.7, phi = 2.1, gamma1 = -6.3)
dat <- test[[1]] 
dat.id <- test[[2]] 

# Get the event rate for the  data set
table(dat.id$event)[2] / (table(dat.id$event)[1] + table(dat.id$event)[2])

# Plot the data
xyplot(y ~ time | event, group=id, data = dat, type = "l")

################################
## Set empty lists for output ##
################################

n.iter = 5000
numbersims = 1 # 100
Event   <- vector("list", numbersims)
Patient <- vector("list", numbersims)
Biomarker <- vector("list", numbersims)

Event3   <- vector("list", numbersims)
Patient3 <- vector("list", numbersims)
Biomarker3 <- vector("list", numbersims)


SumI   <- vector("list", numbersims)
SumII   <- vector("list", numbersims)
SumIII   <- vector("list", numbersims)

AUCI   <- vector("list", numbersims)
AUCII   <- vector("list", numbersims)
AUCIII   <- vector("list", numbersims)

PEI   <- vector("list", numbersims)
PEII   <- vector("list", numbersims)
PEIII   <- vector("list", numbersims)





# Simulate validation data

# Make smaller to speed up process
testv <- Simulate(n = 1000, mean.Cens = 3.75, phi = 2.5, gamma1 = -7.5)
datv <- testv[[1]] 
datv.id <- testv[[2]] 





# Start the clock!
start_time <- Sys.time()
######################
##      SIMULATE    ##
######################

for (count in 1:numbersims) {
  set.seed(13017 + count)
  
  ##############
  ## GET DATA ##
  ##############
  # Simulate development data
  test <- Simulate(n = 2000, mean.Cens = 4.5, phi = 2.5, gamma1 = -9)
  dat <- test[[1]] 
  dat.id <- test[[2]] 

  # Save characteristics of data
  Event[[count]] <- table(dat.id$event)[2] / (table(dat.id$event)[1] + table(dat.id$event)[2])
  Patient[[count]] <- nrow(dat.id)
  Biomarker[[count]] <- nrow(dat)
  
  
  ##############
  ##  MODEL I ##
  ##############

  # Input MODEL I
  lmeObject <- lme(y ~ ns(time, 2) * group, data = dat, 
                   control = lmeControl(opt="optim"),
                   random = list(id = pdDiag(form = ~ ns(time, 2))))
  survObject <- coxph(Surv(Time, event) ~ group , data = dat.id, x = TRUE)
  

  # Obtain elements for JAGS and fit JAGS model
  Data1 <- PrepJAGS(survObject, lmeObject, timeVar = "time")


  # Get model from model text file
  outI <- jags(data = Data1, parameters.to.save = params, model.file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/Code/2017/08-06/model2.txt",
               n.chains = 3, n.adapt = n.iter/2, n.iter = n.iter, n.burnin = n.iter/2,
               n.thin = 2, store.data = TRUE, parallel = TRUE, n.cores = 4)

  # Save the parameters
  SumI[[count]] <- outI$summary[c(1:10),]

  # Make JM object of JAGS
  Extra1 <- ExtractExtra(survObject, lmeObject, timeVar = "time")
  outPI <- JAGSTOJM(outI, Extra1)
  

  # Obtain AUC's
  AUCI[[count]] <-  c(aucJM(outPI, datv, Tstart = 1, Thoriz = 3)$auc, 
                      aucJM(outPI, datv, Tstart = 3, Thoriz = 6)$auc)
  
  PEI[[count]] <- c(prederrJM(outPI, newdata = datv, Tstart = 1,  Thoriz = 3, idVar = "id", interval = FALSE)$prederr,
             prederrJM(outPI, newdata = datv, Tstart = 3,  Thoriz = 6, idVar = "id", interval = FALSE)$prederr)

  ##############
  ## Make CC  ##
  ##############
  
  # Make case-cohort data
  CCdata <- CaseCohort(dat = dat, dat.id = dat.id, event = "event", id = "id", y = "y")
  
  datcc2 <- CCdata[[1]]
  datcc3 <- CCdata[[2]]
  datcc3.id <- CCdata[[3]]
  
  Event3[[count]] <- table(datcc3.id$event)[2] / (table(datcc3.id$event)[1] + table(datcc3.id$event)[2])
  Patient3[[count]] <- nrow(datcc3.id)
  Biomarker3[[count]] <- nrow(datcc3)
  
  ##############
  ## MODEL II ##
  ##############
  
  # Get Data2 by changing the y variable in Data1
  Data2 <- Data1
  Data2$y <- (datcc2$y2)


  # Fit the JAGS model
  outII <- jags(data = Data2, parameters.to.save = params, model.file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/Code/2017/08-06/model2.txt",
                n.chains = 3, n.adapt = n.iter/5, n.iter = n.iter, n.burnin = n.iter/5,
                n.thin = 2, store.data = TRUE, parallel = TRUE, n.cores = 3)


  SumII[[count]] <- outII$summary[c(1:10),]
  
 rm(outII1, outII2, outII3, outII4, outII5)
  # Make JM object of JAGS
  outPII <- JAGSTOJM(outII, Extra1)

    # Obtain AUC's
  AUCII[[count]] <-  c(aucJM(outPII, datv, Tstart = 1, Thoriz = 3)$auc, 
                       aucJM(outPII, datv, Tstart = 3, Thoriz = 6)$auc)
    # Obtain Prediction errors
  PEII[[count]] <- c(prederrJM(outPII, newdata = datv, Tstart = 1,  Thoriz = 3, idVar = "id", interval = FALSE)$prederr,
                       prederrJM(outPII, newdata = datv, Tstart = 3,  Thoriz = 6, idVar = "id", interval = FALSE)$prederr)
  
  ###############
  ## MODEL III ##
  ###############
  
    # Input MODEL III
  lmeObject3 <- lme(y ~ ns(time, 2) * group, data = datcc3, 
                    random = list(id = pdDiag(form = ~ ns(time, 2))))
  survObject3 <- coxph(Surv(Time, event) ~ group , data = datcc3.id, x = TRUE)
  
  
  # Obtain elements for JAGS and fit JAGS model
  Data3 <- PrepJAGS(survObject3, lmeObject3, timeVar = "time")

  # Fit the JAGS model
  outIII <- jags(data = Data3, parameters.to.save = params, model.file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/Code/2017/08-06/model2.txt",
                 n.chains = 3, n.adapt = n.iter/2, n.iter = n.iter, n.burnin = n.iter/2,
                 n.thin = 2, store.data = TRUE, parallel = TRUE, n.cores = 4)
  
  SumIII[[count]] <- outIII$summary[c(1:10),]
  
  # Make JM object of JAGS
  Extra3 <- ExtractExtra(survObject3, lmeObject3, timeVar = "time")
  outPIII <- JAGSTOJM(outIII, Extra3)
 
   # Obtain AUC's
  AUCIII[[count]] <-  c(aucJM(outPIII, datv, Tstart = 1, Thoriz = 3)$auc, 
                        aucJM(outPIII, datv, Tstart = 3, Thoriz = 6)$auc)
  # Obtain Prediction errors
  PEIII[[count]] <- c(prederrJM(outPIII, newdata = datv, Tstart = 1,  Thoriz = 3, idVar = "id", interval = FALSE)$prederr,
                       prederrJM(outPIII, newdata = datv, Tstart = 3,  Thoriz = 6, idVar = "id", interval = FALSE)$prederr)
  
  
}

# Stop the clock
end_time <- Sys.time()
end_time - start_time
#########################################################
## When not all simulations are run give them a number ##
#########################################################

Event1 <- Event
Event31 <- Event3

Patient1 <- Patient
Patient31 <- Patient3

Biomarker1 <- Biomarker
Biomarker31 <- Biomarker3

SumI1 <- SumI
SumII1 <- SumII
SumIII1 <- SumIII

AUCI1 <- AUCI
AUCII1 <- AUCII
AUCIII1 <- AUCIII

PEI1 <- PEI
PEII1 <- PEII
PEIII1 <- PEIII

outI1 <- outI
outII1 <- outII
outIII1 <- outIII

#################
## Save output ##
#################

save(Event1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Event1.RData")
save(Event31,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Event31.RData")

save(Patient1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Patient1.RData")
save(Patient31,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Patient31.RData")

save(Biomarker1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Biomarker1.RData")
save(Biomarker31,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/Biomarker31.RData")

save(SumI1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/SumI1.RData")
save(SumII1,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/SumII1.RData")
save(SumIII1, file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/SumIII1.RData")

save(AUCI1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/AUCI1.RData")
save(AUCII1,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/AUCII1.RData")
save(AUCIII1, file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/AUCIII1.RData")

save(PEI1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/PEI1.RData")
save(PEII1,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/PEII1.RData")
save(PEIII1, file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/PEIII1.RData")

save(outI1,   file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/outI1.RData")
save(outII1,  file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/outII1.RData")
save(outIII1, file = "V:/Users/024017x/Case Cohort/JAGS/Simulations/FINAL/05-03/Output/outIII1.RData")

beep()


# Assess convergence
par(mfrow = c(4,2))
traceplot(outII, parameters = c("betas", "alphas", "gammas"))

