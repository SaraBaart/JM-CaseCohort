
JAGSTOJM <- function(outI, Extra1){
  

# Initialize NULL list
object_JAGS1 <- list(NULL)

# Attach dens Functions
object_JAGS1$Funs$densLong <- function (y, eta.y, scale, log = FALSE, data) { 
  dnorm(x = y, mean = eta.y, sd = scale, log = log)
}

object_JAGS1$Funs$densRE <- function (b, mu, D = NULL, invD = NULL, log = FALSE, prop = TRUE) {
  JMbayes:::dmvnorm(b, mu = mu, Sigma = D, invSigma = invD, log = log, 
                    prop = prop) 
}
#----------------------------------------------------------------------------------------------

# Correct class of the object
class(object_JAGS1) <- "JMbayes"

# ADD Event/Censoring Times (NOTE: "Extra1" is due to the RData object you sent me, specify accordingly)
object_JAGS1$y$Time <- Extra1$Time

# ADD name of timeVar (NOTE: "Extra1" is due to the RData object you sent me, specify accordingly)
object_JAGS1$timeVar <- Extra1$timeVar

# ADD Terms (NOTE: "Extra1" is due to the RData object you sent me, specify accordingly) 
object_JAGS1$Terms$termsYx <- Extra1$TermsX
object_JAGS1$Terms$termsYz <- Extra1$TermsZ
object_JAGS1$Terms$termsT <- Extra1$TermsT
object_JAGS1$Forms$formYz <- formula(Extra1$TermsZ)

# ADD posterior Mean estimates (NOTE: outI refers to the JAGS object)
object_JAGS1$postMeans$betas <- outI$mean$betas
object_JAGS1$postMeans$sigma <- 1/(sqrt(outI$mean$tau))
object_JAGS1$postMeans$D <- solve(outI$mean$inv.D)
object_JAGS1$postMeans$gammas <- outI$mean$gammas
object_JAGS1$postMeans$alphas <- outI$mean$alphas
object_JAGS1$postMeans$Bs.gammas <- outI$mean$Bs.gammas
object_JAGS1$postMeans$Dalphas <- outI$mean$Dalphas
object_JAGS1$postMeans$shapes <- outI$mean$shapes



#---------------------------------------ADD MCMC samples---------------------------------------------------- 
object_JAGS1$mcmc$betas <- as.matrix(data.frame("(Intercept)" = c(outI$samples[[1]][, "betas[1]"], 
                                                                 outI$samples[[2]][, "betas[1]"], 
                                                                 outI$samples[[3]][, "betas[1]"]),  
                                               "ns(year, 2)1" = c(outI$samples[[1]][, "betas[2]"], 
                                                                  outI$samples[[2]][, "betas[2]"], 
                                                                  outI$samples[[3]][, "betas[2]"]) , 
                                               "ns(year, 2)2" = c(outI$samples[[1]][, "betas[3]"], 
                                                                  outI$samples[[2]][, "betas[3]"], 
                                                                  outI$samples[[3]][, "betas[3]"]) , 
                                               "group1" = c(outI$samples[[1]][, "betas[4]"], 
                                                               outI$samples[[2]][, "betas[4]"], 
                                                               outI$samples[[3]][, "betas[4]"]) , 
                                               "ns(year, 2)1:group1" = c(outI$samples[[1]][, "betas[5]"], 
                                                                            outI$samples[[2]][, "betas[5]"], 
                                                                            outI$samples[[3]][, "betas[5]"]) , 
                                               "ns(year, 2)2:group1" = c(outI$samples[[1]][, "betas[6]"], 
                                                                            outI$samples[[2]][, "betas[6]"], 
                                                                            outI$samples[[3]][, "betas[6]"])))

object_JAGS1$mcmc$sigma <- as.matrix(data.frame("sigma" = 1/(sqrt(c(outI$samples[[1]][, "tau"], 
                                                                   outI$samples[[2]][, "tau"], 
                                                                   outI$samples[[3]][, "tau"])))))

object_JAGS1$mcmc$gammas <- as.matrix(data.frame("group" = c(outI$samples[[1]][, "gammas"], 
                                                                outI$samples[[2]][, "gammas"], 
                                                                outI$samples[[3]][, "gammas"])))
object_JAGS1$mcmc$Dalphas <- NULL
object_JAGS1$mcmc$shapes <- NULL
object_JAGS1$mcmc$Bs.gammas <- as.matrix(data.frame("Bs.gammas1" =  c(outI$samples[[1]][, "Bs.gammas[1]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[1]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[1]"]), 
                                                   "Bs.gammas2" =  c(outI$samples[[1]][, "Bs.gammas[2]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[2]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[2]"]), 
                                                   "Bs.gammas3" =  c(outI$samples[[1]][, "Bs.gammas[3]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[3]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[3]"]), 
                                                   "Bs.gammas4" =  c(outI$samples[[1]][, "Bs.gammas[4]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[4]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[4]"]), 
                                                   "Bs.gammas5" =  c(outI$samples[[1]][, "Bs.gammas[5]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[5]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[5]"]), 
                                                   "Bs.gammas6" =  c(outI$samples[[1]][, "Bs.gammas[6]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[6]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[6]"]), 
                                                   "Bs.gammas7" =  c(outI$samples[[1]][, "Bs.gammas[7]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[7]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[7]"]), 
                                                   "Bs.gammas8" =  c(outI$samples[[1]][, "Bs.gammas[8]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[8]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[8]"]), 
                                                   "Bs.gammas9" =  c(outI$samples[[1]][, "Bs.gammas[9]"], 
                                                                     outI$samples[[2]][, "Bs.gammas[9]"], 
                                                                     outI$samples[[3]][, "Bs.gammas[9]"])))



head(object_JAGS1$mcmc$D.inv)
object_JAGS1$mcmc$D.inv <- as.matrix(data.frame("D[1, 1]" = c(outI$samples[[1]][, "inv.D[1,1]"], 
                                                             outI$samples[[2]][, "inv.D[1,1]"], 
                                                             outI$samples[[3]][, "inv.D[1,1]"]), 
                                               "D[2, 1]" = c(outI$samples[[1]][, "inv.D[2,1]"], 
                                                             outI$samples[[2]][, "inv.D[2,1]"], 
                                                             outI$samples[[3]][, "inv.D[2,1]"]) , 
                                               "D[3, 1]" = c(outI$samples[[1]][, "inv.D[3,1]"], 
                                                             outI$samples[[2]][, "inv.D[3,1]"], 
                                                             outI$samples[[3]][, "inv.D[3,1]"]) , 
                                               "D[1, 2]" = c(outI$samples[[1]][, "inv.D[1,2]"], 
                                                             outI$samples[[2]][, "inv.D[1,2]"], 
                                                             outI$samples[[3]][, "inv.D[1,2]"]) , 
                                               "D[2, 2]" = c(outI$samples[[1]][, "inv.D[2,2]"], 
                                                             outI$samples[[2]][, "inv.D[2,2]"], 
                                                             outI$samples[[3]][, "inv.D[2,2]"]) , 
                                               "D[3, 2]" = c(outI$samples[[1]][, "inv.D[3,2]"], 
                                                             outI$samples[[2]][, "inv.D[3,2]"], 
                                                             outI$samples[[3]][, "inv.D[3,2]"]) , 
                                               "D[1, 3]" = c(outI$samples[[1]][, "inv.D[1,3]"], 
                                                             outI$samples[[2]][, "inv.D[1,3]"], 
                                                             outI$samples[[3]][, "inv.D[1,3]"]) , 
                                               "D[2, 1]" = c(outI$samples[[1]][, "inv.D[2,1]"], 
                                                             outI$samples[[2]][, "inv.D[2,1]"], 
                                                             outI$samples[[3]][, "inv.D[2,1]"]) , 
                                               "D[3, 3]" = c(outI$samples[[1]][, "inv.D[3,3]"], 
                                                             outI$samples[[2]][, "inv.D[3,3]"], 
                                                             outI$samples[[3]][, "inv.D[3,3]"])))
# INVERT THE MCMC D VALUES
D.cor <- matrix(rep(NA, times = length(object_JAGS1$mcmc$D.inv)), ncol = ncol(object_JAGS1$mcmc$D.inv), 
                nrow = nrow(object_JAGS1$mcmc$D.inv))
for(i in 1:nrow(object_JAGS1$mcmc$D.inv)){
  solv <- solve(matrix(object_JAGS1$mcmc$D.inv[1, ], ncol = 3, nrow = 3))
  D.cor[i, ] <- as.vector(solv) 
}
#---------------------------------------------------------------------------------------------------------
object_JAGS1$mcmc$alphas <- as.matrix(data.frame("Assoct" = c(outI$samples[[1]][, "alphas"], 
                                                             outI$samples[[2]][, "alphas"], 
                                                             outI$samples[[3]][, "alphas"])))
#---------------------------------------------FINISH MCMC SAMPLES------------------------------------------------------------


# these take default values if not specified otherwise
object_JAGS1$control$GQsurv <- "GaussKronrod"
object_JAGS1$control$GQsurv.k <- outI$data$K
object_JAGS1$y$lag <- outI$data$lag
object_JAGS1$control$ordSpline <- 4
object_JAGS1$param <- "td-value"
object_JAGS1$estimateWeightFun <- FALSE
object_JAGS1$Funs$transFun.value <- function(x, data) {
  x
}
object_JAGS1$Funs$hasScale <- TRUE

# PLacement of knots for baseline risk function (usually returned by JMbayesObject but not from JAGSObject
# so we input manually) 
object_JAGS1$control$knots <- outI$data$kn

# Longitudinal observed measurements
object_JAGS1$y$y <- outI$data$y

#rm(D.cor, solv, Extra1, i, outI)
 return(object_JAGS1)

}
