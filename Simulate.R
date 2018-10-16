Simulate <- function(n) {


# Try to simulate the data based on this
n <- n # number of subjects 
K <- 25  # number of planned repeated measurements per subject, per outcome
t.max <- 15 # maximum follow-up time 

################################################

# parameters for the linear mixed effects model
betas <- c("Group0" = 1, "Group1" = 0.1, "Time1" = 0.3,
           "Time1:Time1" = 0.1)

sigma.y <- 1 # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" = -7.5, "Group" = -2) # coefficients for baseline covariates 
alpha <- 1 # association parameter 
phi <- 2 # shape for the Weibull baseline hazard 
mean.Cens <- 3.2 # mean of the exponential distribution for the censoring mechanism 

D <- diag(c(0.99, 0.99, 0.07)^2)


################################################

Bkn <- c(0.5, t.max + 0.5)
kn <- c(mean.Cens)


# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
DF <- data.frame(obstime = times, gender = factor(rep(group, each = K)))
X <- model.matrix(~ gender + obstime + I(obstime*obstime), data = DF)
Z <- model.matrix(~ obstime + I(obstime*obstime), data = DF) 

# design matrix for the survival model
W <- cbind("(Intercept)" = 1, "Group" = group)


################################################

#simulate random effects
b <- mvrnorm(n, rep(0, nrow(D)), D)


# simulate longitudinal responses
id <- rep(1:n, each = K)
eta.y <- as.vector(X %*% betas + rowSums(Z * b[id, ])) # linear predictor
y <- rnorm(n * K, eta.y, sigma.y)

# simulate event times
eta.t <- as.vector(W %*% gammas)
invS <- function (t, u, i) {
  h <- function (s) {
    group0 <- 1 - group[i]
    group1 <- group[i]
    XX <- cbind(group0, group1, s, (s^2))
    ZZ <- cbind(1, s, s^2)
    f1 <- as.vector(XX %*% betas + rowSums(ZZ * b[rep(i, nrow(ZZ)), ]))
    exp(log(phi) + (phi - 1) * log(s) + eta.t[i] + f1 * alpha)
  }
  integrate(h, lower = 0, upper = t)$value + log(u)
}
u <- runif(n)
trueTimes <- numeric(n)
for (i in 1:n) {
  Up <- 50
  tries <- 20
  Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  while(inherits(Root, "try-error") && tries > 0) {
    tries <- tries - 1
    Up <- Up + 200
    Root <- try(uniroot(invS, interval = c(1e-05, Up), u = u[i], i = i)$root, TRUE)
  }
  trueTimes[i] <- if (!inherits(Root, "try-error")) Root else NA
}
na.ind <- !is.na(trueTimes)
trueTimes <- trueTimes[na.ind]
W <- W[na.ind, , drop = FALSE]
long.na.ind <- rep(na.ind, each = K)
y <- y[long.na.ind]
X <- X[long.na.ind, , drop = FALSE]
Z <- Z[long.na.ind, , drop = FALSE]
DF <- DF[long.na.ind, ]
n <- length(trueTimes)

# simulate censoring times from an exponential distribution,
# and calculate the observed event times, i.e., min(true event times, censoring times)
Ctimes <- runif(n, 0, 2 * mean.Cens)
Time <- pmin(trueTimes, Ctimes)
event <- as.numeric(trueTimes <= Ctimes) # event indicator


################################################


# keep the nonmissing cases, i.e., drop the longitudinal measurements
# that were taken after the observed event time for each subject.
ind <- times[long.na.ind] <= rep(Time, each = K)
y <- y[ind]
X <- X[ind, , drop = FALSE]
Z <- Z[ind, , drop = FALSE]
id <- id[long.na.ind][ind]
id <- match(id, unique(id))
id.id <- unique(id)

dat <- DF[ind, ]
dat$id <- id
dat$y <- y
dat$Time <- Time[id]
dat$event <- event[id]
dat.id <- data.frame(id = id.id, Time = Time, event = event, group = W[,2])
names(dat) <- c("time", "group", "id", "y", "Time", "event")



# Imitate a case-cohort design
sample.ids <- sample(unique(dat$id), replace = FALSE, size = round(0.3333*length(unique(dat$id))))
event.ids <-  dat$id[ which(dat$event ==1)]
cc.ids <- unique(c(sample.ids, event.ids) ) 
dat$y2 <- ifelse(dat$id %in% cc.ids, dat$y, NA)   
dat$CC <- ifelse(dat$id %in% cc.ids, 1,0) 
           
           
# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W,
   betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, mean.Cens, t.max,
   trueTimes, u, Root, invS, D, b, K,
   times, group, i, tries, Up, Bkn, kn, DF, id.id)

return(dat)

}
