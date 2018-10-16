
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

listbetas <- paste0("betas[", c(1:ncol(outI$data$X)),"]")
test <- vector("list", length(listbetas))
  for (i in 1:length(listbetas)){
test[[i]]  <- data.frame(x = c(outI$samples[[1]][, listbetas[[i]] ],
                               outI$samples[[2]][, listbetas[[i]] ],
                               outI$samples[[3]][, listbetas[[i]] ]))
}
   object_JAGS1$mcmc$betas <- as.matrix(as.data.frame(test))
   colnames(object_JAGS1$mcmc$betas) <-  colnames(outI$data$X)


object_JAGS1$mcmc$sigma <- as.matrix(data.frame("sigma" = 1/(sqrt(c(outI$samples[[1]][, "tau"], 
                                                                   outI$samples[[2]][, "tau"], 
                                                                   outI$samples[[3]][, "tau"])))))

listgammas <- paste0("gammas[", c(1:ncol(outI$data$W)),"]")
test <- vector("list", length(listgammas))
for (i in 1:length(listgammas)){
  test[[i]]  <- data.frame(x = c(outI$samples[[1]][, listgammas[[i]] ],
                                 outI$samples[[2]][, listgammas[[i]] ],
                                 outI$samples[[3]][, listgammas[[i]] ]))
}
object_JAGS1$mcmc$gammas <- as.matrix(as.data.frame(test))
colnames(object_JAGS1$mcmc$gammas) <-  colnames(outI$data$W)


object_JAGS1$mcmc$Dalphas <- NULL
object_JAGS1$mcmc$shapes <- NULL

ncol(outI$data$W2)
listBsgammas <- paste0("Bs.gammas[", c(1:ncol(outI$data$W2)),"]")
test <- vector("list", length(listBsgammas))
for (i in 1:length(listBsgammas)){
  test[[i]]  <- data.frame(x = c(outI$samples[[1]][, listBsgammas[[i]] ],
                                 outI$samples[[2]][, listBsgammas[[i]] ],
                                 outI$samples[[3]][, listBsgammas[[i]] ]))
}
object_JAGS1$mcmc$Bs.gammas <- as.matrix(as.data.frame(test))
colnames(object_JAGS1$mcmc$Bs.gammas) <-  paste0("Bs.gammas", c(1:ncol(outI$data$W2)))




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
