# 30th March 2022
# This file runs the three real life examples of:
# Kerr - Irish Rath n = 255 obs
# Buchanan - Palaeo-Indian demography n = 628 obs
# Armit - Population Decline in Iron Age n=2021 obs
# using the new IntCal20 curve

# Read in IntCal20 curve
calcurve <- read.table("Curves/intcal20.14c", sep = ",", header=FALSE, skip=11)
names(calcurve) <- c("calage", "c14age", "c14sig", "Delta14C", "DeltaSigma")

# Read in the necessary functions
source('WalkerDirichletMixtureUpdateFunsFinal.R') # This also reads in the slice sampling SliceUpdateFuns.R
source('NealDirichletMixtureMasterFunctionsFinal.R')
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")

# Choose number of iterations for sampler
niter <- 100000 
nthin <- 10 # Ensure that (niter/nthin)/2 is not too small (this is max posterior sample) 

##### Read in the data
ExampleSet <- "Armit"

if(ExampleSet == "Kerr") {
  Kerr <- read.csv("RealDatasets/kerr2014sss_sup.csv", header = FALSE, sep =  ",")
  x <- Kerr[,3]
  xsig <- Kerr[,4]
} else if(ExampleSet == "Buchanan") {
  Data <- read.csv("RealDatasets/buchanan2008pde.csv", header = FALSE, sep =",")
  x <- Data[,3]
  xsig <- Data[,4]
} else if(ExampleSet == "Armit") {
  Armit <- read.csv("RealDatasets/armit2014rcc_sd01.csv", header = TRUE, sep = ",")
  # Remove the two observations with missing xsig values
  remove <- which(is.na(Armit$error))
  Armit <- Armit[-remove, ]
  x <- Armit$X14C.age
  xsig <- Armit$error  
} else {
  stop("Unknown Example")
}

############################################################################
# Now choose hyperparameters
############################################################
# Prior on the concentration parameter
# Place  a gamma prior on alpha
# alpha ~ Gamma(alphaprshape, alphaprrate) 
# A small alpha means more concentrated (i.e. few clusters)
# Large alpha not concentrated (many clusters)
cprshape <- alphaprshape <- 1 
cprrate <- alphaprrate <- 1 

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

# Run the Walekr slice sampler 
WalkerTemp <- WalkerBivarDirichlet(x = x, xsig = xsig, 
                                   lambda = lambda, nu1 = nu1, nu2 = nu2, 
                                   A = A, B = B, 
                                   cprshape = cprshape, cprrate = cprrate, 
                                   niter = niter, nthin = nthin, theta = inittheta,
                                   slicew = max(1000, diff(range(x))/2), m = 10, calcurve = calcurve, kstar = 10)

##############################
# Also find the SPD estimate to plot alongside
##############################
# Find the independent calibration probabilities
yrange <- floor(range(WalkerTemp$theta))
yfromto <- seq(max(0,yrange[1]-400), min(50000, yrange[2]+400), by = 1)

# Find the calibration curve mean and sd over the yrange
CurveR <- FindCalCurve(yfromto, calcurve)

# Now we want to apply to each radiocarbon determination
# Matrix where each column represents the posterior probability of each theta in yfromto 
indprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = CurveR$curvemean, calsig = CurveR$curvesd))

# Find the SPD estimate (save as dataframe)
SPD <- data.frame(calage = yfromto,
                   prob = apply(indprobs, 1, sum)/dim(indprobs)[2])

source("RealExamplesWalkerPostProcessingFinal.R")

if(ExampleSet == "Kerr") {
  save.image(file = "KerrNPBayesWalker.RData")
} else if(ExampleSet == "Buchanan") {
  save.image(file = "BuchananNPBayes.RData")
} else if(ExampleSet == "Armit") {
  save.image(file = "ArmitNPBayes.RData")
}


