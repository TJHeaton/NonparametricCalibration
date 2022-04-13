# 5th April 2022
# Uniform Distribution Example:
# Range restriced to Holocene
# Edited so that the choice of lambda, nu1 and nu2 are adaptive to inittheta
# Corrected Walker Sampler 

# In this file we will try a simulation study to show how assuming a shared underlying density improves calibration
# The simulation study will proceed by:
# Assuming a true underlying uniform density
# Drawing N observations from those clusters and calculating corresponding radiocarbon determinations 
# Calibrating these N observations in two scenarios:
# a - Entirely independently
# b - Assuming a shared underlying density

# In this version we use the updated (simplified) JRSSB code for the samplers 
# This is different from the initial submission since the Walker version is corrected

set.seed(81) # For revised paper


# Read in the necessary functions
source('WalkerDirichletMixtureUpdateFunsFinal.R') # This also reads in the slice sampling SliceUpdateFuns.R
source('NealDirichletMixtureMasterFunctionsFinal.R')
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")


# Read in IntCal20 curve
calcurve <- read.table("Curves/intcal20.14c", sep = ",", header=FALSE, skip=11)
names(calcurve) <- c("calage", "c14age", "c14sig", "Delta14C", "DeltaSigma")


# Choose number of iterations for sampler
niter <- 10000 # 10000 # 10000
nthin <- 5

# Create limits for the Uniform distribution
startrange <- c(100,14000)
phaserange <- c(50, 1000)
nclus <- 1

############################################################################
# Now choose fixed DP hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(alphaprshape, alphaprrate) 
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- alphaprshape <- 1 
cprrate <- alphaprrate <- 1 

detno <- c(50, 100, 200, 500)
nrep <- 50  # 50 # 20

Storex1 <- matrix(NA, nrow = nrep, ncol = 50)
Storetheta1 <- matrix(NA, nrow = nrep, ncol = 50)
Storephasestart1 <- matrix(NA, nrow = nrep, ncol = nclus)
Storephaselength1 <- matrix(NA, nrow = nrep, ncol = nclus)


Storex2 <- matrix(NA, nrow = nrep, ncol = 100)
Storetheta2 <- matrix(NA, nrow = nrep, ncol = 100)
Storephasestart2 <- matrix(NA, nrow = nrep, ncol = nclus)
Storephaselength2 <- matrix(NA, nrow = nrep, ncol = nclus)


Storex3 <- matrix(NA, nrow = nrep, ncol = 200)
Storetheta3 <- matrix(NA, nrow = nrep, ncol = 200)
Storephasestart3 <- matrix(NA, nrow = nrep, ncol = nclus)
Storephaselength3 <- matrix(NA, nrow = nrep, ncol = nclus)


Storex4 <- matrix(NA, nrow = nrep, ncol = 500)
Storetheta4 <- matrix(NA, nrow = nrep, ncol = 500)
Storephasestart4 <- matrix(NA, nrow = nrep, ncol = nclus)
Storephaselength4 <- matrix(NA, nrow = nrep, ncol = nclus)


WalkerNPLoss <- array(NA, dim = c(nrep, 2, length(detno)), dimnames = list(1:nrep, c("l1loss","l2loss"), detno))
NealNPLoss <- array(NA, dim = c(nrep, 2, length(detno)), dimnames = list(1:nrep, c("l1loss","l2loss"), detno))
IndCalLoss <- array(NA, dim = c(nrep, 2, length(detno)), dimnames = list(1:nrep, c("l1loss","l2loss"), detno))

# Do the longest runs n = 500 first 
for(k in length(detno):1) {
  nobs <- detno[k]
  # Determine sd on the radiocarbon determinations
  xsig <- rep(25, nobs) # Assume an sd of 25 on the radiocarbon determination
  cat("nobs = ", nobs, "\n")
  for(i in 1:nrep) {
    # Create simulated data restricted to (100,49500)
    inrange <- FALSE
    while(!inrange) { # Edited to allow different weights for mixture distribution 
      SimData <- CreateSingleUniformPhaseData(nobs, startrange = startrange, calcurve, phaserange = phaserange, xsig)
      inrange <- (max(SimData$thetatrue) < 49500) & (min(SimData$thetatrue) > 100) 
      cat(".")
    }  
    thetatrue <-SimData$thetatrue
    x <- SimData$x
    if(k == 1) {
      Storex1[i,] <- x
      Storetheta1[i,] <- thetatrue
      Storephasestart1[i,] <- SimData$phasestart
      Storephaselength1[i,] <- SimData$phaselength
    }
    if(k == 2) {
      Storex2[i,] <- x
      Storetheta2[i,] <- thetatrue
      Storephasestart2[i,] <- SimData$phasestart
      Storephaselength2[i,] <- SimData$phaselength
    }
    if(k == 3) {
      Storex3[i,] <- x
      Storetheta3[i,] <- thetatrue
      Storephasestart3[i,] <- SimData$phasestart
      Storephaselength3[i,] <- SimData$phaselength
    }
    if(k == 4) {
      Storex4[i,] <- x
      Storetheta4[i,] <- thetatrue
      Storephasestart4[i,] <- SimData$phasestart
      Storephaselength4[i,] <- SimData$phaselength
    }
    
    #### Updated adaptive version
    # Prior on mu theta for DP - very uninformative based on observed data
    initprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = calcurve$c14age, calsig = calcurve$c14sig))
    inittheta <- calcurve$calage[apply(initprobs, 2, which.max)]
    # Choose A and B from range of theta
    A <- median(inittheta)
    B <- 1 / (max(inittheta) - min(inittheta))^2
    maxrange <- max(inittheta) - min(inittheta)
    
    # Parameters for sigma2 (sigma^2 ~ InvGamma(nu1, nu2))
    # E[tau] = (1/100)^2 Var[tau] = (1/100)^4
    # Interval for sigma2 is approx 1/ c(nu2/nu1 - 2*nu2^2/nu1, nu2/nu1 + 2*nu2^2/nu1)
    tempspread <- 0.1 * mad(inittheta)
    tempprec <- 1/(tempspread)^2
    nu1 <- 0.25
    nu2 <- nu1 / tempprec
    
    # Setup the NP method
    lambda <- (100/maxrange)^2  # Each muclust ~ N(mutheta, sigma2/lambda)
    
    #####################################################################
    # Implement the Neal version of the DPMM
    NealTemp <- BivarGibbsDirichletwithSlice(x = x, xsig = xsig,
                                             lambda = lambda, nu1 = nu1, nu2 = nu2,
                                             A = A, B = B,
                                             mualpha = NA, sigalpha = NA,
                                             alphaprshape = alphaprshape, alphaprrate = alphaprrate,
                                             niter = niter, nthin = nthin, theta = inittheta,
                                             w = max(1000, diff(range(x))/2), m = 10,
                                             calcurve = calcurve, nclusinit = 10, showprogress = FALSE)
    
    NealNPLoss[i,,k] <- unlist(FindLoss(NealTemp$theta, true = thetatrue))
    cat("+")
    
    #####################################################################
    # Implement the Walker version of the DPMM
    WalkerTemp <- WalkerBivarDirichlet(x = x, xsig = xsig, 
                                       lambda = lambda, nu1 = nu1, nu2 = nu2, 
                                       A = A, B = B, 
                                       cprshape = cprshape, cprrate = cprrate, 
                                       niter = niter, nthin = nthin, theta = inittheta,
                                       slicew = max(1000, diff(range(x))/2), m = 10, 
                                       calcurve = calcurve, kstar = 10, showprogress = FALSE)
    
    WalkerNPLoss[i,,k] <- unlist(FindLoss(WalkerTemp$theta, true = thetatrue))
    cat("*")
    
    #####################################################################
    #### Now perform independent calibration of the determinations
    yrange <- floor(range(WalkerTemp$theta))
    yfromto <- seq( max(0,yrange[1]-400), min(50000, yrange[2]+400), by = 1)
    
    # Find the calibration curve mean and sd over the yrange
    CurveR <- FindCalCurve(yfromto, calcurve)
    
    # Now we want to apply to each radiocarbon determination
    # Matrix where each column represents the posterior probability of each theta in yfromto 
    indprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = CurveR$curvemean, calsig = CurveR$curvesd))
    
    # Now apply function which works out the mean posterior loss
    IndLoss <- sapply(1:length(x), function(i, Mprobs, Vtheta, yfromto) 
      FindIndLoss(Mprobs[,i], true = Vtheta[i], yfromto = yfromto), Mprobs = indprobs, Vtheta = thetatrue, yfromto = yfromto)
    
    # Now find the mean loss over all determinations
    IndCalLoss[i,,k] <- apply(IndLoss, 1, mean)
  }
  cat("***\n")
}


outdigits <- 1
# Print the output into a LaTeX table
library(xtable)
NPLoss <- NealNPLoss
Improvement <- 100 * (1 - (NPLoss/IndCalLoss))
SimSummaryMean <- t(apply(Improvement, c(2,3), mean)) # This gives the mean reduction in loss by using NP Bayes (it is pretty good)
SimSummaryImprove <- t(apply(Improvement, c(2,3), function(x) mean(x > 0)))
SimSummaryImproveCI <- t(apply(Improvement, c(2,3), function(x) 
  paste("(", paste(100 * round(binom.test(sum(x > 0), length(x))$conf.int, digits = 2), collapse = ","), ")", sep = "")))

SimSummaryBest <- t(apply(Improvement, c(2,3), max))
SimSummaryWorst <- t(apply(Improvement, c(2,3), min))

cols.combined <- 5 * ncol(SimSummaryMean)
SimCombinedOut <-  matrix(NA, nrow=nrow(SimSummaryMean), ncol = cols.combined )
SimCombinedOut[, seq(1, cols.combined, 5)] <- as.integer(100 * SimSummaryImprove)
SimCombinedOut[, seq(2, cols.combined, 5)] <- SimSummaryImproveCI
SimCombinedOut[, seq(3, cols.combined, 5)] <- round(SimSummaryMean, digits = outdigits)
SimCombinedOut[, seq(4, cols.combined, 5)] <- round(SimSummaryBest, digits = outdigits)
SimCombinedOut[, seq(5, cols.combined, 5)] <- round(SimSummaryWorst, digits = outdigits)
colnames(SimCombinedOut) <- rep(c("Prop. improved (%)", "CI", "Mean", "Maximum", "Minimum"), 2)

rownames(SimCombinedOut) <- detno

NealSimCombinedOut <- SimCombinedOut

#######################################################################

Improvement <- 100 * (1 - (WalkerNPLoss/IndCalLoss))
SimSummaryMean <- t(apply(Improvement, c(2,3), mean)) # This gives the mean reduction in loss by using NP Bayes (it is pretty good)
SimSummaryImprove <- t(apply(Improvement, c(2,3), function(x) mean(x > 0)))
SimSummaryImproveCI <- t(apply(Improvement, c(2,3), function(x) 
  paste("(", paste(100 * round(binom.test(sum(x > 0), length(x))$conf.int, digits = 2), collapse = ","), ")", sep = "")))

SimSummaryBest <- t(apply(Improvement, c(2,3), max))
SimSummaryWorst <- t(apply(Improvement, c(2,3), min))

cols.combined <- 5 * ncol(SimSummaryMean)
SimCombinedOut <-  matrix(NA, nrow=nrow(SimSummaryMean), ncol = cols.combined )
SimCombinedOut[, seq(1, cols.combined, 5)] <- as.integer(100 * SimSummaryImprove)
SimCombinedOut[, seq(2, cols.combined, 5)] <- SimSummaryImproveCI
SimCombinedOut[, seq(3, cols.combined, 5)] <- round(SimSummaryMean, digits = outdigits)
SimCombinedOut[, seq(4, cols.combined, 5)] <- round(SimSummaryBest, digits = outdigits)
SimCombinedOut[, seq(5, cols.combined, 5)] <- round(SimSummaryWorst, digits = outdigits)
colnames(SimCombinedOut) <- rep(c("Prop. improved (%)", "CI", "Mean", "Maximum", "Minimum"), 2)

rownames(SimCombinedOut) <- detno

WalkerSimCombinedOut <- SimCombinedOut

BothDPMMSimOut <- cbind(NealSimCombinedOut, WalkerSimCombinedOut)


xtable(BothDPMMSimOut)
#digits(xtmp)[4:7] <- 1
#print(xtmp, include.rownames = FALSE)

# # Find a CI for the proportion where there is an improvement
# apply(Improvement, c(2,3), function(x) paste(round(binom.test(sum(x > 0), length(x))$conf.int, digits = 2)))
# ImproveCI <- t(apply(Improvement, c(2,3), function(x) 
#   paste("(", paste(round(binom.test(sum(x > 0), length(x))$conf.int, digits = 2), collapse = ","), ")", sep = "")))
# #Produce output in a table
# xtable(SimCombinedOut, digits = c(1,rep(c(0, rep(3, 3)), 2)))

save.image("SimOutput/WalkerandNealSimStudyUniformRestrictedRange_Revisedv1AdaptiveSeed81.RData")

