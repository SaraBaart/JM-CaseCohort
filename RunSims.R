####################
## Load libraries ## 
####################
library(JMbayes)
library(jagsUI)
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
setwd("V:/Users/024017x/Case Cohort/Simulaties/05-03")
# Simulate data
source("Functions/Simulate.R")
# Prepare Data for output
source("Functions/PrepJAGS.R")
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
# Simulate data
test <- Simulate(n = 2000, mean.Cens = 2.5, phi = 2, sigma.y = 1,
                 gamma1 = -9.5,
                 gamma2 = -2,
                 beta1 = 1,
                 beta2 = 0.1,
                 beta3 = 0.3,
                 beta4 = 0.1,
                 alpha = 1)
dat <- test[[1]] 
dat.id <- test[[2]] 
nrow(dat.id)/2000
# Get the event rate for the  data set
table(dat.id$event)[2] / (table(dat.id$event)[1] + table(dat.id$event)[2])

# Plot the data
xyplot(y ~ time | event, group=id, data = dat, type = "l")




################################
## Set empty lists for output ##
################################

n.iter = 7000
numbersims = 200 # 100
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








######################
##      SIMULATE    ##
######################
start_time <- Sys.time()
for (count in 1:numbersims) {
  set.seed(20181101 + count)
  
  ##############
  ## GET DATA ##
  ##############
  # Simulate development data
  test <- Simulate(n = 2000, mean.Cens = 2.5, phi = 2, sigma.y = 1,
                   gamma1 = -9.5,
                   gamma2 = -2,
                   beta1 = 1,
                   beta2 = 0.1,
                   beta3 = 0.3,
                   beta4 = 0.1,
                   alpha = 1)
  dat <- test[[1]] 
  dat.id <- test[[2]] 

  
  # Simulate validation data
  testv <-  Simulate(n = 2000, mean.Cens = 2.5, phi = 2, sigma.y = 1,
                     gamma1 = -9.5,
                     gamma2 = -2,
                     beta1 = 1,
                     beta2 = 0.1,
                     beta3 = 0.3,
                     beta4 = 0.1,
                     alpha = 1)
  datv <- testv[[1]] 
  datv.id <- testv[[2]] 
  
  
  # Save characteristics of data
  Event[[count]] <- table(dat.id$event)[2] / 
    (table(dat.id$event)[1] + table(dat.id$event)[2])
  Patient[[count]] <- nrow(dat.id)
  Biomarker[[count]] <- nrow(dat)
  
  
  ##############
  ##  MODEL I ##
  ##############

  # Input MODEL I
   lmeObject <- lme(y ~ group + time + I(time^2), data = dat, 
                   control = lmeControl(opt="optim"),
                   random = list(id = pdDiag(form = ~ time + I(time^2))))
  
  survObject <- coxph(Surv(Time, event) ~ group , data = dat.id, x = TRUE)

  
  # Obtain elements for JAGS and fit JAGS model
  Data1 <- PrepJAGS(survObject, lmeObject, timeVar = "time")


  # Get model from model text file
  outI <- jags(data = Data1, parameters.to.save = params, 
                   model.file = "Model/model.txt", n.chains = 3, 
                   n.adapt = n.iter/5, n.iter = n.iter, n.burnin = n.iter/5,
                   n.thin = 5, store.data = TRUE, parallel = TRUE, n.cores = 4)

  # Save the parameters
  SumI[[count]] <- outI$summary[c(1:10),]

  # Make JM object of JAGS
  Extra1 <- ExtractExtra(survObject, lmeObject, timeVar = "time")
  outPI <- JAGSTOJM(outI, Extra1)
  

  # Obtain AUC's
  AUCI[[count]] <-  c(aucJM(outPI, datv, Tstart = 1, Thoriz = 2)$auc, 
                      aucJM(outPI, datv, Tstart = 1.5, Thoriz = 3)$auc)
  
  PEI[[count]] <- c(prederrJM(outPI, newdata = datv, Tstart = 1,  
                    Thoriz = 2, idVar = "id", interval = FALSE)$prederr,
             prederrJM(outPI, newdata = datv, Tstart = 1.5,  
                    Thoriz = 3, idVar = "id", interval = FALSE)$prederr)

  ##############
  ## Make CC  ##
  ##############
  
  # Make case-cohort data
  CCdata <- CaseCohort(dat = dat, dat.id = dat.id, event = "event", 
                       id = "id", y = "y")
  
  datcc2 <- CCdata[[1]]
  datcc3 <- CCdata[[2]]
  datcc3.id <- CCdata[[3]]
  
  Event3[[count]] <- table(datcc3.id$event)[2] / 
    (table(datcc3.id$event)[1] + table(datcc3.id$event)[2])
  Patient3[[count]] <- nrow(datcc3.id)
  Biomarker3[[count]] <- nrow(datcc3)
  
  ##############
  ## MODEL II ##
  ##############
  
  # Get Data2 by changing the y variable in Data1
  Data2 <- Data1
  Data2$y <- (datcc2$y2)

  # Fit the JAGS model

  outII <- jags(data = Data2, parameters.to.save = params, 
                model.file = "Model/model.txt", n.chains = 3, 
                n.adapt = n.iter, n.iter = n.iter*2, n.burnin = n.iter,
                n.thin = 5, store.data = TRUE, parallel = TRUE, n.cores = 4)
  

  SumII[[count]] <- outII$summary[c(1:10),]
  
 
  # Make JM object of JAGS
  outPII <- JAGSTOJM(outII, Extra1)

    # Obtain AUC's
  AUCII[[count]] <-  c(aucJM(outPII, datv, Tstart = 1, Thoriz = 2)$auc, 
                       aucJM(outPII, datv, Tstart = 1.5, Thoriz = 3)$auc)
    # Obtain Prediction errors
  PEII[[count]] <-  c(prederrJM(outPII, newdata = datv, Tstart = 1,  
                      Thoriz = 2, idVar = "id", interval = FALSE)$prederr,
                      prederrJM(outPII, newdata = datv, Tstart = 1.5,  
                      Thoriz = 3, idVar = "id", interval = FALSE)$prederr)
  
  ###############
  ## MODEL III ##
  ###############
  
    # Input MODEL III
  lmeObject3 <- lme(y ~ group + time + I(time^2), data = datcc3, 
                    random = list(id = pdDiag(form = ~ time + I(time^2))))
  survObject3 <- coxph(Surv(Time, event) ~ group , data = datcc3.id, x = TRUE)
  

  
  # Obtain elements for JAGS and fit JAGS model
  Data3 <- PrepJAGS(survObject3, lmeObject3, timeVar = "time")

  # Fit the JAGS model
  outIII <- jags(data = Data3, parameters.to.save = params, 
                 model.file = "Model/model.txt", n.chains = 3, 
                 n.adapt = n.iter/5, n.iter = n.iter, n.burnin = n.iter/5,
                 n.thin = 5, store.data = TRUE, parallel = TRUE, n.cores = 4)
  
  
 
  
  SumIII[[count]] <- outIII$summary[c(1:10),]
  
  # Make JM object of JAGS
  Extra3 <- ExtractExtra(survObject3, lmeObject3, timeVar = "time")
  outPIII <- JAGSTOJM(outIII, Extra3)
 
   # Obtain AUC's
  AUCIII[[count]] <-   c(aucJM(outPIII, datv, Tstart = 1, Thoriz = 2)$auc, 
                        aucJM(outPIII, datv, Tstart = 1.5, Thoriz = 3)$auc)
  # Obtain Prediction errors
  PEIII[[count]] <-  c(prederrJM(outPIII, newdata = datv, Tstart = 1,  
                       Thoriz = 2, idVar = "id", interval = FALSE)$prederr,
                       prederrJM(outPIII, newdata = datv, Tstart = 1.5,  
                       Thoriz = 3, idVar = "id", interval = FALSE)$prederr)
  
  
}
end_time <- Sys.time()
end_time - start_time
beep()


#########################################################
## When not all simulations are run give them a number ##
#########################################################
Event
Event2 <- Event
Event23 <- Event3

Patient2 <- Patient
Patient23 <- Patient3

Biomarker2 <- Biomarker
Biomarker23 <- Biomarker3

SumI2 <- SumI
SumII2 <- SumII
SumIII2 <- SumIII

AUCI2 <- AUCI
AUCII2 <- AUCII
AUCIII2 <- AUCIII

PEI2 <- PEI
PEII2 <- PEII
PEIII2 <- PEIII

outI2 <- outI
outII2 <- outII
outIII2 <- outIII

#################
## Save output ##
#################
getwd()
save(Event,   file = "Output/Event.RData")
save(Event3,  file = "Output/Event3.RData")

save(Patient,   file = "Output/Patient.RData")
save(Patient3,  file = "Output/Patient3.RData")

save(Biomarker,   file = "Output/Biomarker.RData")
save(Biomarker3,  file = "Output/Biomarker3.RData")

save(SumI,   file = "Output/SumI.RData")
save(SumII,  file = "Output/SumII.RData")
save(SumIII, file = "Output/SumIII.RData")

save(AUCI,   file = "Output/AUCI.RData")
save(AUCII,  file = "Output/AUCII.RData")
save(AUCIII, file = "Output/AUCIII.RData")

save(PEI,   file = "Output/PEI.RData")
save(PEII,  file = "Output/PEII.RData")
save(PEIII, file = "Output/PEIII.RData")

save(outI,   file = "Output/outI.RData")
save(outII,  file = "Output/outII.RData")
save(outIII, file = "Output/outIII.RData")


# Assess convergence
par(mfrow = c(2,2))
traceplot(outIII, parameters = c("betas", "alphas", "gammas", "inv.D"))

str(PEI)
str(PEI1)
 SumI[[61]][5,]
 SumI1[[61]][5,]
outI$summary[1:7,]
outII$summary[1:7,]
outIII$summary[1:7,]

 