# This code will plot the predictive density for Walker

Temp <- WalkerTemp

out <- Temp

# Find predictive density
npost <- dim(Temp$delta)[1]
nburn <- floor(npost/2)
npostsum <- 5000

# Create the range over which to plot the density 
# Note since tempx is not seq(T_a, T_b, 1) (i.e. spaced every year) then
# sum(postdenCI) != 1 
# We should have sum(postdenCI) * gridgap = 1 (where gridgap is the spacing in tempx)
tempx <- seq(floor(min(Temp$theta, na.rm = TRUE)), ceiling(max(Temp$theta, na.rm = TRUE)), by = 1) 

#Now choose the ids of the posterior sample
sampid <- sample(x = nburn:npost, size = npostsum, replace = npostsum > (npost-nburn))

# Now work out the actual posterior predictive density
postDmat <- apply(as.row(sampid), 2, function(i, out, x, lambda, nu1, nu2) 
  WalkerFindpred(x, w = out$w[[i]], phi = out$phi[[i]], tau = out$tau[[i]], 
                 muphi = out$muphi[i], lambda = lambda, nu1 = nu1, nu2 = nu2), 
  out = Temp, x = tempx, lambda = lambda, nu1 = nu1, nu2 = nu2)

# Find CI and mean
postdenCI <- apply(postDmat, 1, quantile, probs = c(0.025, 0.975))
postden <- apply(postDmat, 1, mean)

# Create a pdf output file according to the example
oname <- paste("OutputRealExamples/", ExampleSet, "NPDensityWalkerVersionAdaptive.pdf", sep = "") 

# Open a pdf file for the plot 
scaleplot <- 0.8
pdf(oname, width = scaleplot * 12, height = scaleplot * 8)
# pdf(oname, width = scaleplot * 12, height = scaleplot * 6) # Old vesion - too small to see density estimate

# Create a layout which has 2/3 of the plot showing the predictive density and 1/3 showing the number of clusters
layout.matrix <- matrix(c(1, 2), nrow = 1, ncol = 2)
layout(mat = layout.matrix,
       heights = c(1), # Heights of the two rows
       widths = c(10, 4.5)) # Widths of the two columns

SPDcol <- grey(0.1, alpha = 0.5)

if(ExampleSet == "Kerr") {
  calcurve$ADages <- 1950.5 - calcurve$calage
  ADtempx <- 1950.5 - tempx 
  SPD$ADages <- 1950.5 - SPD$calage
  
  # Plot the predictive joint density 
  par(mar = c(5, 4.5, 2, 2) + 0.1, las = 1)
  plot(calcurve$ADages, calcurve$c14age, col = "blue", 
       ylim = range(x) + c(-2,2) * max(xsig), xlim = range(ADtempx),
       xlab = "Calendar Year (AD)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
       type = "l")
  calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
  calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
  lines(calcurve$ADages, calcurve$ub, lty = 2, col = "blue" )
  lines(calcurve$ADages, calcurve$lb, lty = 2, col = "blue")
  polygon(c(rev(calcurve$ADages), calcurve$ADages), c(rev(calcurve$lb), calcurve$ub), 
          col=rgb(0,0,1,.3), border=NA)
  rug(x, side = 2)
  
  
  par(new = TRUE, las = 1)
  plot(ADtempx, postden, lty = 1, col = "purple", 
       ylim = c(0, 2*max(postdenCI)), xlim = range(ADtempx),
       axes = FALSE, xlab = NA, ylab = NA, type = "n")
  polygon(c(SPD$ADages, rev(SPD$ADages)), c(SPD$prob, rep(0, length(SPD$prob))), 
          border = NA, col = SPDcol)
  lines(ADtempx, postden, col = "purple")
  lines(ADtempx, postdenCI[1,], col = "purple", lty = 2)
  lines(ADtempx, postdenCI[2,], col = "purple", lty = 2)
  #axis(side = 4, at = axTicks(4)/2, col = "purple", col.ticks = "purple")
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
        line = -1.1)

  legend("topright", legend = c("IntCal20", "NP Density Estimate - Walker", "SPD Estimate" ), 
         lty = c(1,1,-1), pch = c(NA, NA, 15), 
         col = c("blue", "purple", SPDcol), cex = 0.8)


# Plot the number of clusters
  WalkerNClust <- apply(WalkerTemp$delta, 1, function(x) length(unique(x)))
  WalkerNClust <- WalkerNClust[nburn:npost]
  hist(WalkerNClust, xlab = "Number of Clusters", main = "",
       probability = TRUE, breaks = seq(0.5, max(WalkerNClust)+1, by = 1))
  mtext(paste0("(", letters[2], ")"), side = 3, adj = 1, 
        line = -1.1)
} else { # Just do a calBP plot for Buchanan and Armit
  if(ExampleSet == "Buchanan") {
    callimits <- c(16000, 8000) 
    ylims <- c(5000, 13500)
    xaxs <- "i"
    denscale <- 1.15
  } else {
    callimits <- rev(range(tempx))
    ylims <- range(x) + c(-2,2) * quantile(xsig, 0.9)
    xaxs<- "i"
    denscale <- 1.4
  }
    
  # Plot the predictive joint density 
  par(mar = c(5, 4.5, 2, 2) + 0.1, las = 1)
  plot(calcurve$calage, calcurve$c14age, col = "blue", 
       ylim = ylims, xlim = callimits,
       xlab = "Calendar Age (cal yr BP)", ylab = expression(paste(""^14, "C", " age (", ""^14, "C yr BP)")),
       type = "l", xaxs = xaxs)
  calcurve$ub <- calcurve$c14age + 1.96 * calcurve$c14sig
  calcurve$lb <- calcurve$c14age - 1.96 * calcurve$c14sig
  lines(calcurve$calage, calcurve$ub, lty = 2, col = "blue" )
  lines(calcurve$calage, calcurve$lb, lty = 2, col = "blue")
  polygon(c(rev(calcurve$calage), calcurve$calage), c(rev(calcurve$lb), calcurve$ub), 
          col=rgb(0,0,1,.3), border=NA)
  rug(x, side = 2)
  
  
  par(new = TRUE, las = 1)
  plot(tempx, postden, lty = 1, col = "purple",
       ylim = c(0, denscale*max(postdenCI)), xlim = callimits,
       axes = FALSE, xlab = NA, ylab = NA, type = "n")
  polygon(c(SPD$calage, rev(SPD$calage)), c(SPD$prob, rep(0, length(SPD$prob))), 
          border = NA, col = SPDcol)
  lines(tempx, postden, col = "purple")
  lines(tempx, postdenCI[1,], col = "purple", lty = 2)
  lines(tempx, postdenCI[2,], col = "purple", lty = 2)
  #axis(side = 4, at = axTicks(4)/2, col = "purple", col.ticks = "purple")
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
        line = -1.1)
  
  # Add a polygon for the Glen West Bog in Armit
  if(ExampleSet == "Armit") {
    Bogx <- c(1950.5 + 786, 1950.5 + 703)
    polygon(c(Bogx, rev(Bogx)), c(-0.01,-0.01, 3*max(postdenCI),3*max(postdenCI)), 
            col = rgb(1, 0, 0, .2), border = NA)
    # This is setup so that red shaded box will cover whole plot
  }
  if(ExampleSet == "Buchanan") {
    YDx <- c(12807, 11653)
    polygon(c(YDx, rev(YDx)), c(-0.01,-0.01, 3*max(postdenCI),3*max(postdenCI)), 
            col = rgb(1, 0, 0, .2), border = NA)
    text(x = mean(YDx), y = max(postdenCI), labels = "YD", adj = 0.5)
  }
  
  mtext(paste0("(", letters[1], ")"), side = 3, adj = 0.05, 
        line = -1.1)
  
  
  legend("topright", legend = c("IntCal20", "NP Density Estimate - Walker", "SPD Estimate" ), 
         lty = c(1,1,-1), pch = c(NA, NA, 15), 
         col = c("blue", "purple", SPDcol), cex = 0.8)
  
  # Plot the number of clusters - in plot (b)
  WalkerNClust <- apply(WalkerTemp$delta, 1, function(x) length(unique(x)))
  WalkerNClust <- WalkerNClust[nburn:npost]
  hist(WalkerNClust, xlab = "Number of Clusters", main = "",
       probability = TRUE, breaks = seq(0.5, max(WalkerNClust)+1, by = 1))
  mtext(paste0("(", letters[2], ")"), side = 3, adj = 1, 
        line = -1.1)
}

# Now close the plot
dev.off()



