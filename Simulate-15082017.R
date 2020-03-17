Simulate <- function(n, mean.Cens, phi, gamma1) {
 

# Try to simulate the data based on this
n <- n # number of subjects 1000
K <- 25  # number of planned repeated measurements per subject, per outcome
t.max <- 15 # maximum follow-up time 14.10579

################################################

# parameters for the linear mixed effects model
betas <- c("Group0" = 1, "Group1" = 0.1, "Group0:Time1" = 1.2,
           "Group1:Time1" = 1, "Group0:Time2" = 1.3, "Group1:Time2" = -0.5)

sigma.y <- 0.3 # measurement error standard deviation


# parameters for the survival model
gammas <- c("(Intercept)" = gamma1, "Group" = 1) # coefficients for baseline covariates 0.1324
alpha <- 1 #0.3 # association parameter 1.1818
phi <- phi # shape for the Weibull baseline hazard 2.284594 1.884594
mean.Cens <- mean.Cens # mean of the exponential distribution for the censoring mechanism 4.5

D <- diag(c(0.7, 0.7, 0.7)^2)


################################################

Bkn <- c(0, 10.10579)
kn <- c(5.031541)


# design matrices for the longitudinal measurement model
# but this can be easily generalized
times <- c(replicate(n, c(0, sort(runif(K-1, 0, t.max))))) # at which time points longitudinal measurements are supposed to be taken
group <- rep(0:1, each = n/2) # group indicator, i.e., '0' placebo, '1' active treatment
DF <- data.frame(obstime = times, gender = factor(rep(group, each = K)))
X <- model.matrix(~ gender*ns(obstime, knots = kn, Boundary.knots = Bkn), data = DF)
Z <- model.matrix(~ ns(obstime, knots = kn, Boundary.knots = Bkn), data = DF)

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
    NS <- ns(s, knots = kn, Boundary.knots = Bkn)
    XX <- cbind(group0, group1, group0*NS[, 1], group1*NS[, 1],
                group0*NS[, 2], group1*NS[, 2])
    ZZ <- cbind(1, NS)
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

#summary(tapply(id, id, length))
#table(event)
#n
#mean(event)

# Make a spaghetti plot of the data
xyplot(y ~ time | event, group=id, data=dat, type = "l")

# Count the number of events
table(dat.id$event)[2]/ (table(dat.id$event)[1] + table(dat.id$event)[2])


# delete all unused objects
rm(y, X, Z, id, n, na.ind, long.na.ind, ind, Ctimes, Time, event, W,
   betas, sigma.y, gammas, alpha, eta.t, eta.y, phi, mean.Cens, t.max,
   trueTimes, u, Root, invS, D, b, K,
   times, group, i, tries, Up, Bkn, kn, DF, id.id)

return(list(dat, dat.id))

}
