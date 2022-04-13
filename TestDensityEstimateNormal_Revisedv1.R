# 5th April 2022

# Study to see how well the densities are reconstructed
set.seed(14)


# Read in the necessary functions
source('WalkerDirichletMixtureUpdateFunsFinal.R') # This also reads in the slice sampling SliceUpdateFuns.R
source('NealDirichletMixtureMasterFunctionsFinal.R')
source("WalkerMasterFunctionFinal.R")
source("SimStudyFuncsFinal.R")


# Read in IntCal20 curve
calcurve <- read.table("Curves/intcal20.14c", sep = ",", header=FALSE, skip=11)
names(calcurve) <- c("calage", "c14age", "c14sig", "Delta14C", "DeltaSigma")

# Number of observations and xsig
nobs <- 100
xsig <- rep(25, nobs)

# Choose base f(theta)
wtrue <- c(0.1, 0.4, 0.5)
wtrue <- wtrue / sum (wtrue)
phitrue <- c(3500, 4200, 5000)
tautrue <- 1/c(200,100,300)^2
nclus <- length(wtrue)

# Create some observed data from the clusters according to probabilities
ctrue <- sample(1:nclus, nobs, replace = TRUE, prob = wtrue)
thetatrue <- rnorm(nobs, mean = phitrue[ctrue], sd = 1/sqrt(tautrue[ctrue]))
hist(thetatrue, breaks = 20)

#### Now create some radiocarbon determinations x

# Interpolate calobration curve mean and sd at theta values
calinterp <- FindCal(thetatrue, calmu = calcurve$c14age, caltheta = calcurve$calage, calsig = calcurve$c14sig)

# Sample some calibration curve values
xcalcurve <- rnorm(nobs, calinterp$mu, calinterp$sigma)

x <- rnorm(nobs, mean = xcalcurve, sd = xsig)  

###############################################
# Now run the Walker Sampler
###############################################
# Choose number of iterations for sampler
niter <- 50000 # 10000 # 10000
nthin <- 5

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
                                         calcurve = calcurve, nclusinit = 10)

#####################################################################
# Implement the Walker version of the DPMM
WalkerTemp <- WalkerBivarDirichlet(x = x, xsig = xsig, 
                                   lambda = lambda, nu1 = nu1, nu2 = nu2, 
                                   A = A, B = B, 
                                   cprshape = cprshape, cprrate = cprrate, 
                                   niter = niter, nthin = nthin, theta = inittheta,
                                   slicew = max(1000, diff(range(x))/2), m = 10, 
                                   calcurve = calcurve, kstar = 10)
Temp <- WalkerTemp

##############################
# Find the independent calibration probabilities
yrange <- floor(range(WalkerTemp$theta))
yfromto <- seq( max(0,yrange[1]-400), min(50000, yrange[2]+400), by = 1)

# Find the calibration curve mean and sd over the yrange
CurveR <- FindCalCurve(yfromto, calcurve)

# Now we want to apply to each radiocarbon determination
# Matrix where each column represents the posterior probability of each theta in yfromto 
indprobs <- mapply(calibind, x, xsig, MoreArgs = list(calmu = CurveR$curvemean, calsig = CurveR$curvesd))


# Find predictive density
npost <- dim(Temp$delta)[1]
nburn <- floor(npost/2)
npostsum <- 5000

# Create the range over which to plot the density
tempx <- seq(floor(min(Temp$theta, na.rm = TRUE)), ceiling(max(Temp$theta, na.rm = TRUE)), by = 1) 

#Now choose the ids of the posterior sample
sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost-nburn))

#############################
# Now work out the actual posterior predictive density
#############################
# Walker method - slice sampling
WalkerpostDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) 
  WalkerFindpred(x, w = out$w[[i]], phi = out$phi[[i]], tau = out$tau[[i]], 
                 muphi = out$muphi[i], lambda = lambda, nu1 = nu1, nu2 = nu2), 
  out = WalkerTemp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2)

WalkerpostdenCI <- apply(WalkerpostDmat, 1, quantile, probs = c(0.025, 0.975))
Walkerpostden <- apply(WalkerpostDmat, 1, mean)

# Neal method - Polya Urn
# Create a matrix where each column is the density for a particular sample id
# We can then find the mean along each row 
NealpostDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) 
  NealFindpred(x, c = out$c[i,], phi = out$phi[[i]], tau = out$tau[[i]], 
           alpha = out$alpha[i], muphi = out$muphi[i], 
           lambda = lambda, nu1 = nu1, nu2 = nu2), 
  out = NealTemp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2)

NealpostdenCI <- apply(NealpostDmat, 1, quantile, probs = c(0.025, 0.975))
Nealpostden <- apply(NealpostDmat, 1, mean)





# Find and plot true exact density
truedensv2 <- function(t, w, truemean, trueprec) {
  dens <- 0
  for(i in 1:length(w)) {
    dens <- dens + w[i] * dnorm(t, mean = truemean[i], sd = 1/sqrt(trueprec[i]))
  }
  dens
}
curve(truedensv2(x, w = wtrue, truemean = phitrue, trueprec = tautrue), 
      from = min(thetatrue) - 200, to = max(thetatrue) + 200, n = 401, 
      ylim = c(0, 0.003), 
      xlim = rev(c(min(thetatrue) - 200, max(thetatrue) + 200)),
      col = "red")


# hist(thetatrue, prob = TRUE, add = TRUE)
lines(rev(tempx), rev(Walkerpostden), col = "purple")
lines(tempx, WalkerpostdenCI[1,], col = "purple", lty = 3)
lines(tempx, WalkerpostdenCI[2,], col = "purple", lty = 3)  

lines(rev(tempx), rev(Nealpostden), col = "blue")
lines(tempx, NealpostdenCI[1,], col = "blue", lty = 3)
lines(tempx, NealpostdenCI[2,], col = "blue", lty = 3) 


pdf("TestPlotNormalDensityRecon.pdf", width = 10, height = 6)

layout.matrix <- matrix(c(1, 1, 2, 3), nrow = 2, ncol = 2)

layout(mat = layout.matrix,
       heights = c(3, 3), # Heights of the two rows
       widths = c(10, 4.5)) # Widths of the two columns

par(mar = c(5, 4.5, 2, 2) + 0.1, las = 1)
ywid <- max(x) - min(x)
plot(calcurve$calage, calcurve$c14age, col = "blue", 
     ylim = range(x) + c(-0.5*ywid,200), rev(c(min(thetatrue) - 200, max(thetatrue) + 200)),
     xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
     type = "l")
calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
lines(calcurve$calage, calcurve$ub, lty = 3, col = "blue" )
lines(calcurve$calage, calcurve$lb, lty = 3, col = "blue")
polygon(c(rev(calcurve$calage), calcurve$calage), c(rev(calcurve$lb), calcurve$ub), col=rgb(0,0,1,.3), border=NA)
rug(x, side = 2)

# Now plot as a rug along the bottom the calendar age estimates
par(new = TRUE, las = 0)
curve(truedensv2(x, w = wtrue, truemean = phitrue, trueprec = tautrue), 
      from = min(thetatrue) - 200, to = max(thetatrue) + 200, n = 401, 
      ylim = 2 * c(0, 0.0024), # So takes up 1/3 of height
      xlim = rev(c(min(thetatrue) - 200, max(thetatrue) + 200)),
      col = "red", axes = FALSE, xlab = NA, ylab = NA, lwd = 2)
lines(rev(tempx), rev(Walkerpostden), col = "purple")
lines(tempx, WalkerpostdenCI[1,], col = "purple", lty = 2)
lines(tempx, WalkerpostdenCI[2,], col = "purple", lty = 2)  

lines(rev(tempx), rev(Nealpostden), col = "blue")
lines(tempx, NealpostdenCI[1,], col = "blue", lty = 2)
lines(tempx, NealpostdenCI[2,], col = "blue", lty = 2) 

# Overlay the independent SPD
SPD <- apply(indprobs, 1, sum)/dim(indprobs)[2]
lines(yfromto, SPD, col = "orange")

legend("topright", legend = c(expression(paste("True density ",f(theta))), "Walker DP", "Neal DP", "Independent SPD"), 
       lty = 1, col = c("red", "purple", "blue", "orange"))
mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
      line = -1.1)

# Also look at the number of clusters
NealNClust <- apply(NealTemp$c[nburn:npost,], 1, max)
WalkerNClust <- apply(WalkerTemp$delta, 1, function(x) length(unique(x)))
WalkerNClust <- WalkerNClust[nburn:npost]
maxclust <- max(c(NealNClust, WalkerNClust))+1

hist(NealNClust, xlab = "Number of Clusters", main = "Neal - P\u{F2}lya Urn DP", 
     probability = TRUE, ylim = c(0, 0.4), 
     breaks = seq(1.5, maxclust, by = 1))
mtext(paste0("(", letters[2], ")"), side = 3, adj = 1, 
      line = -1.1)

# Also look at the number of clusters

hist(WalkerNClust, xlab = "Number of Clusters", main = "Walker - Slice Sample DP",
     probability = TRUE, ylim = c(0, 0.4), 
     breaks = seq(1.5, maxclust, by = 1))
mtext(paste0("(", letters[3], ")"), side = 3, adj = 1, 
      line = -1.1)

dev.off()

