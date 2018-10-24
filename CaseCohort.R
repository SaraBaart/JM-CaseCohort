##########################
##    Load libraries    ## 
##########################

library(JMbayes)
library(jagsUI)
library(splines)
library(MASS)
library(nlme)

####################################
##    Simulate and prepare data   ##
####################################


# Simulate data with missing biomarker values (y2) for patients outside the case cohort design, using the function Simulate.
data <- Simulate(n = 2000)

# For the longitudinal submodel; plug in the general mean of y2 in the new longitudinal variable y3.
data$y3 <- ifelse(is.na(data$y2), mean(data$y2, na.rm = TRUE), data$y2)

# Make a short version of the data to fit the survival submodel
data.id  <- data[!duplicated(data$id),]

# Fit the longitudinal submodel
lmeObject <- lme(y3 ~ group + time + I(time^2), data = data, 
                 random = list(id = pdDiag(form = ~ time + I(time^2))))

# Fit the survival submodel
survObject <- coxph(Surv(Time, event) ~ group , data = data.id, x = TRUE)


# Obtain elements needed to fit the joint model with JAGS. 
# The timeVar is the time variable in the longitudinal submodel.
Data <- PrepJAGS(survObject, lmeObject, timeVar = "time")

# Replace the longitudinal variable y with the version with missings (y2)
Data$y <- (data$y2)


################################
##       Write the model      ##
################################


# Write the joint model in BUGS code
writeLines("
 model{
    for (i in 1:N) {
    # Longitudinal Part
    for (j in offset[i]:(offset[i+1] - 1)) {
      muy[j] <- inprod(b[i, 1:ncZ], Z[j, 1:ncZ])
      
      y[j] ~ dnorm(muy[j], tau)
    }
    
    # Survival Part
    etaBaseline[i] <- inprod(gammas[1], W[i])
    log.h0.T[i] <-  inprod(Bs.gammas[1:ncW2], W2[i, 1:ncW2])
    f.T[i] <- inprod(b[i, 1:ncZ], Ztime[i, 1:ncZ])
    log.hazard[i] <- log.h0.T[i] + etaBaseline[i] + alphas * f.T[i]
    for (k in 1:K) {
      log.h0.s[i, k] <- inprod(Bs.gammas[1:ncW2], W2s[K * (i - 1) + k, 1:ncW2])
      f.s[i, k] <- inprod(b[i, 1:ncZ], Zs[K*(i - 1) + k, 1:ncZ])
      SurvLong[i, k] <- wk[k] * exp(log.h0.s[i, k] + alphas * f.s[i, k])
    }
    log.survival[i] <- - exp(etaBaseline[i]) * P[i] * sum(SurvLong[i, ])
    lambda[i] <- C - (event[i] * log.hazard[i]) - log.survival[i]
    zeros[i] ~ dpois(lambda[i])
    # Random Effects Part
    
    mu_b[i, 1] <- betas[1] + betas[4] * W[i]
    mu_b[i, 2] <- betas[2] 
    mu_b[i, 3] <- betas[3] 
    b[i, 1:nb] ~ dmnorm(mu_b[i, ], inv.D[, ])
  }
  #  Priors
  # Longitudinal Part
  betas[1:ncX] ~ dmnorm(priorMean.betas[], priorTau.betas[, ])
  tau ~ dgamma(priorA.tau, priorB.tau)
  # Survival Part
  gammas[1:ncW] ~ dmnorm(priorMean.gammas[], priorTau.gammas[ ])
  alphas ~ dnorm(priorMean.alphas, priorTau.alphas)
  Bs.gammas[1:ncW2] ~ dmnorm(priorMean.Bs.gammas[], priorTau.Bs.gammas[, ])
  # Random Effects Part
  inv.D[1:nb, 1:nb] ~ dwish(priorR.D[, ], priorK.D)
}
", con = "model.txt")


# Set which parameters to extract
params <- c("betas", "alphas", "gammas", "inv.D", "b", "Bs.gammas", "tau")

#############################
##      Fit the model      ##
#############################

# Fit the model using jags
out <- jags(data = Data, parameters.to.save = params, 
              model.file = "model.txt", n.chains = 3, 
              n.adapt = 1000, n.iter = 5000, n.burnin = 2000,
              n.thin = 5, store.data = TRUE, parallel = TRUE, n.cores = 4)

# Rewrite the output as jointmodelBayes object to calculate predictive accuracy measures
Extra1 <- ExtractExtra(survObject, lmeObject, timeVar = "time")
outPI <- JAGSTOJM(outI, Extra1)
