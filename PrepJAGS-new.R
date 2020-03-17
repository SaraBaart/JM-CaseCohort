  
PrepJAGS <- function(survObject, lmeObject, timeVar){

timeVar <- "time"
  timeVar <- timeVar
  survMod <- "spline-PH"
  param <- "td-value"
  robust = FALSE; robust.b = FALSE; df = 4; df.b = 4; lag = 0
  init = NULL; extraForm = NULL; priors = NULL

  # extract response & design matrix survival process
  TermsT <- survObject$terms
  formT <- formula(survObject)
  W <- survObject$x
  Time <- survObject$y[, 1]
  nT <- length(Time)
  if (!length(W))
    W <- NULL
  event <- survObject$y[, 2]
  zeros <- numeric(nT) # for the zeros trick
  # longitudinal process
  id <- as.vector(unclass(lmeObject$groups[[1]]))
  b <- data.matrix(ranef(lmeObject))
  one.RE <- ncol(b) == 1
  if (one.RE)
    b <- cbind(b, rep(0, nrow(b)))
  dimnames(b) <- NULL
  nY <- nrow(b)
  if (nY != nT)
    stop("sample sizes in the longitudinal and event processes differ.\n")
  TermsX <- lmeObject$terms
  data <- lmeObject$data[all.vars(TermsX)]
  data <- data[complete.cases(data), ]
  formYx <- formula(lmeObject)
  mfX <- model.frame(TermsX, data = data)
  X <- model.matrix(formYx, mfX)
  formYz <- formula(lmeObject$modelStruct$reStruct[[1]])
  mfZ <- model.frame(terms(formYz), data = data)
  TermsZ <- attr(mfZ, "terms")
  Z <- model.matrix(formYz, mfZ)
  if (one.RE)
    Z <- cbind(Z, rep(0, nrow(Z)))
  y.long <- model.response(mfX, "numeric")

  offset <- as.vector(c(1, 1 + cumsum(tapply(id, id, length))))
  
  

  # control values
  con <- list(program = "JAGS", n.chains = 1, n.iter = 10000, n.burnin = 5000,
              n.thin = 1, n.adapt = 1000, K = 100, C = 5000, working.directory = getwd(),
              bugs.directory = "C:/Program Files/WinBUGS14/", openbugs.directory = NULL,
              clearWD = TRUE, over.relax = TRUE, knots = NULL, ObsTimes.knots = TRUE,
              lng.in.kn = 5, ordSpline = 4, bugs.seed = 1, quiet = FALSE)
  ##########
  data.id <- data[!duplicated(id), ]
  if (!timeVar %in% names(data))
    stop("\n'timeVar' does not correspond to one of the columns in the model.frame of 'lmeObject'.")
  # extra design matrices for the longitudinal part
  data.id[[timeVar]] <- pmax(Time - lag, 0)
  if (param %in% c("td-value", "td-both")) {
    mfX.id <- model.frame(TermsX, data = data.id)
    mfZ.id <- model.frame(TermsZ, data = data.id)
    Xtime <- model.matrix(formYx, mfX.id)
    Ztime <- model.matrix(formYz, mfZ.id)
    if (one.RE)
      Ztime <- cbind(Ztime, rep(0, nrow(Ztime)))
    long <- c(X %*% fixef(lmeObject)) + rowSums(Z * b[id, ])
  }
  if (param %in% c("td-extra", "td-both")) {
    mfX.deriv <- model.frame(terms(extraForm$fixed), data = data)
    TermsX.deriv <- attr(mfX.deriv, "terms")
    mfZ.deriv <- model.frame(terms(extraForm$random), data = data)
    TermsZ.deriv <- attr(mfZ.deriv, "terms")
    mfX.deriv.id <- model.frame(TermsX.deriv, data = data.id)
    mfZ.deriv.id <- model.frame(TermsZ.deriv, data = data.id)
    Xtime.deriv <- model.matrix(extraForm$fixed, mfX.deriv.id)
    Ztime.deriv <- model.matrix(extraForm$random, mfZ.deriv.id)
    Xderiv <- model.matrix(extraForm$fixed, mfX.deriv)
    Zderiv <- model.matrix(extraForm$random, mfZ.deriv)
    long.deriv <- as.vector(c(Xderiv %*% fixef(lmeObject)[extraForm$indFixed]) +
                              if (length(extraForm$indRandom) > 1 || extraForm$indRandom)
                                rowSums(Zderiv * b[id, extraForm$indRandom, drop = FALSE])
                            else
                              rep(0, nrow(Zderiv)))
  }
  if (param == "td-value")
    long.deriv <- NULL
  if (param == "td-extra")
    long <- NULL

  # response vectors and design matrices
  y <- list(y = y.long, offset = offset, logT = log(Time), event = event, zeros = zeros, lag = lag)
  x <- list(X = X, Z = Z, W = W)
  x <- switch(param,
              "td-value" = c(x, list(Xtime = Xtime, Ztime = Ztime)),
              "td-extra" = c(x, list(Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
              "td-both" = c(x, list(Xtime = Xtime, Ztime = Ztime,
                                    Xtime.deriv = Xtime.deriv, Ztime.deriv = Ztime.deriv)),
              "shared-RE" = x)
  wk <- JM:::gaussKronrod()$wk
  sk <- JM:::gaussKronrod()$sk
  K <- length(sk)
  P <- Time/2
  st <- outer(P, sk + 1)
  id.GK <- rep(seq_along(Time), each = K)
  data.id2 <- data.id[id.GK, ]
  data.id2[[timeVar]] <- c(t(st))
  x <- c(x, list(P = P, st = st, wk = wk))
  if (param %in% c("td-value", "td-both")) {
    mfX <- model.frame(TermsX, data = data.id2)
    mfZ <- model.frame(TermsZ, data = data.id2)
    Xs <- model.matrix(formYx, mfX)
    Zs <- model.matrix(formYz, mfZ)
    if (one.RE)
      Zs <- cbind(Zs, rep(0, nrow(Zs)))
    x <- c(x, list(Xs = Xs, Zs = Zs))
  }
  if (param %in% c("td-extra", "td-both")) {
    mfX.deriv <- model.frame(TermsX.deriv, data = data.id2)
    mfZ.deriv <- model.frame(TermsZ.deriv, data = data.id2)
    Xs.deriv <- model.matrix(extraForm$fixed, mfX.deriv)
    Zs.deriv <- model.matrix(extraForm$random, mfZ.deriv)
    x <- c(x, list(Xs.deriv = Xs.deriv, Zs.deriv = Zs.deriv))
  }
  evTime <- cens <- Time
  cens[event == 1] <- 0
  evTime[event == 0] <- NA
  y <- c(y, list(evTime = evTime, censTime = cens))
  # extra design matrices for the log approximated baseline hazard
  if (survMod == "spline-PH") {
    kn <- if (is.null(con$knots)) {
      pp <- seq(0, 1, length.out = con$lng.in.kn + 2)
      pp <- tail(head(pp, -1), -1)
      tt <- if (con$ObsTimes.knots) Time else Time[event == 1]
      quantile(tt, pp, names = FALSE)
    } else {
      con$knots
    }
    kn <- kn[kn < max(Time)]
    rr <- sort(c(rep(range(Time, st), con$ordSpline), kn))
    con$knots <- rr
    W2 <- splineDesign(rr, Time, ord = con$ordSpline)
    if (any(colSums(W2) == 0))
      stop("\nsome of the knots of the B-splines basis are set outside the range",
           "\n   of the observed event times for one of the strata; refit the model",
           "\n   setting the control argument 'equal.strata.knots' to FALSE.")
    W2s <- splineDesign(rr, c(t(st)), ord = con$ordSpline)
    x <- c(x, list(W2 = W2, W2s = W2s))
  }
  # default priors
  ncX <- ncol(X)
  ncZ <- ncol(Z)
  ncW <- ncol(x$W)
  ncW2 <- ncol(x$W2)
  betas <- rep(0, ncX)
  var.betas <- rep(con$K, ncX)
   sigma2 <- 1
 # sigma2 <- lmeObject$sigma^2
#  VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
 #              lmeObject$sigma^2)[[1]]
  #VC <- lapply(pdMatrix(lmeObject$modelStruct$reStruct), "*",
   #            sigma2)[[1]]

  # if (one.RE)

  
  #  VC <- cbind(c(VC, 0), c(0, 1))
   VC <- diag(10, nrow = ncZ, ncol = ncZ)
  inv.VC <- solve(VC)
  alphas <- Dalphas <- 0
  var.alphas <- var.Dalphas <- con$K
  gammas <- rep(0, ncW)
  var.gammas <- rep(con$K, ncW)
  Bs.gammas <- rep(0, ncW2)
  var.Bs.gammas <- rep(con$K/10, ncW2)
  # Data to be passed to WinBUGS
  Data <- list(N = nY, df = df, df.b = df.b, C = con$C, K = K,  ncX = ncX,
               ncZ = ncZ, ncW = ncW, ncW2 = ncW2,
               priorMean.betas = betas, priorTau.betas = diag(1 / var.betas),
               priorA.tau = (1/sigma2)^2 / 10, priorB.tau = (1/sigma2) / 10,
               priorMean.gammas = gammas, priorTau.gammas = 1 / var.gammas,
               priorMean.Bs.gammas = Bs.gammas, priorTau.Bs.gammas = diag(1/ var.Bs.gammas),
               priorMean.alphas = alphas, priorTau.alphas = 1 / var.alphas,
               priorMean.Dalphas = alphas, priorTau.Dalphas = 1 / var.Dalphas,
               priorA.sigma.t = 10, priorB.sigma.t = 10, priorA.rho = 2, priorB.rho = 1,
               nb = ncZ, mu0 = rep(0, ncZ), priorR.D = ncZ * inv.VC, priorK.D = ncZ, knots = con$knots)
  Data1 <- c(Data, x, y)

 return(Data1)
}

