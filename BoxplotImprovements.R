# This file takes the output of the simulation study and creates boxplots 
# You must run the simualtion study first as it reads in the .RData files with the output


# Attempt to show boxplot of improvements
library(latex2exp)
# Read in the output from the Three Cluster Simulation Study from where it is saved
# In my case, this is a sibdirectory called SimOutput
load("SimOutput/WalkerandNealSimStudyThreeClusterRestrictedRange_Revisedv1AdaptiveSeed19.RData")


# Find the Walker Improvements 
Walkerl1 <- as_tibble(Improvement[,1,])
Walkerl1$loss <- "l1"

Walkerl2 <- as_tibble(Improvement[,2,])
Walkerl2$loss <- "l2"
WalkerComb <- rbind(Walkerl1, Walkerl2)

WalkerLong <- pivot_longer(WalkerComb,
                        cols = -"loss",
                        names_to = "Size",
                        values_to = "Improvement")
WalkerLong$type <- "Walker Slice DPMM"

# Find the Neal Improvements
Improvement <- 100 * (1 - (NPLoss/IndCalLoss))
Neall1 <- as_tibble(Improvement[,1,])
Neall1$loss <- "l1"

Neall2 <- as_tibble(Improvement[,2,])
Neall2$loss <- "l2"
NealComb <- rbind(Neall1, Neall2)

NealLong <- pivot_longer(NealComb,
                           cols = -"loss",
                           names_to = "Size",
                           values_to = "Improvement")
NealLong$type <- "Neal P\u{F2}lya Urn"

BothMethodsLong <- rbind(WalkerLong, NealLong)

loss.labs <- c(TeX("$l_1$-loss"), TeX("$l_2$-loss"))
loss.labs <- c("l[1]-loss", "l[2]-loss")
names(loss.labs) <- c("l1", "l2")

scalefac <- 1.6
pdf("ThreeClusterSimStudy.pdf", width = scalefac * 33/4, height = scalefac * 0.25 * 47/4)
p3 <- ggplot(BothMethodsLong, aes(x=Size, y=Improvement, fill=type)) +
  geom_boxplot() + facet_wrap(~ loss, labeller = labeller(loss = as_labeller(loss.labs, label_parsed)))
p3 <- p3 + scale_fill_brewer(palette="RdBu") + scale_x_discrete(limits=c("50", "100", "200", "500"))
p3 <- p3 + labs(x = TeX("Sample Size $n$"),
                y = "Percentage Improvement \n (Non-Parametric over Independent Calibration)",
                title = "Example 2: Three underlying normal phases",
                subtitle = TeX("$\\mu_\\phi = 3000$, $\\theta_i \\in \\left[100, 15000\\right]$"))

p3 <- p3 + theme(legend.position = "bottom", 
                 plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
p3 <- p3 + labs(tag = "(b)")
p3 <- p3 + geom_hline(yintercept = 0, col = "blue")
print(p3)
dev.off()

#########################################################
#### Single Cluster Example
rm(list = ls())
load("SimOutput/WalkerandNealSimStudySingleCluster_Revisedv1AdaptiveSeed43.RData")
# Find the Walker Improvements 
Walkerl1 <- as_tibble(Improvement[,1,])
Walkerl1$loss <- "l1"

Walkerl2 <- as_tibble(Improvement[,2,])
Walkerl2$loss <- "l2"
WalkerComb <- rbind(Walkerl1, Walkerl2)

WalkerLong <- pivot_longer(WalkerComb,
                           cols = -"loss",
                           names_to = "Size",
                           values_to = "Improvement")
WalkerLong$type <- "Walker Slice DPMM"

# Find the Neal Improvements
Improvement <- 100 * (1 - (NPLoss/IndCalLoss))
Neall1 <- as_tibble(Improvement[,1,])
Neall1$loss <- "l1"

Neall2 <- as_tibble(Improvement[,2,])
Neall2$loss <- "l2"
NealComb <- rbind(Neall1, Neall2)

NealLong <- pivot_longer(NealComb,
                         cols = -"loss",
                         names_to = "Size",
                         values_to = "Improvement")
NealLong$type <- "Neal P\u{F2}lya Urn"

BothMethodsLong <- rbind(WalkerLong, NealLong)

loss.labs <- c(TeX("$l_1$-loss"), TeX("$l_2$-loss"))
loss.labs <- c("l[1]-loss", "l[2]-loss")
names(loss.labs) <- c("l1", "l2")

scalefac <- 1.6
pdf("SingleClusterSimStudy.pdf", width = scalefac * 33/4, height = scalefac * 0.25 * 47/4)
p1 <- ggplot(BothMethodsLong, aes(x=Size, y=Improvement, fill=type)) +
  geom_boxplot() + facet_wrap(~ loss, labeller = labeller(loss = as_labeller(loss.labs, label_parsed)))
p1 <- p1 + scale_fill_brewer(palette="RdBu") + scale_x_discrete(limits=c("50", "100", "200", "500"))
p1 <- p1 + labs(x = TeX("Sample Size $n$"),
                y = "Percentage Improvement \n (Non-Parametric over Independent Calibration)",
                title = "Example 1: Single underlying normal phase",
                subtitle = TeX("$\\mu_\\phi = 10000$, $\\theta_i \\in \\left[100, 49500\\right]$"))

p1 <- p1 + theme(legend.position = "bottom", 
                 plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))

p1 <- p1 + labs(tag = "(a)")
p1 <- p1 + geom_hline(yintercept = 0, col = "blue")
print(p1)
dev.off()

##########################################################
#### Uniform Example
##########################################################

load("SimOutput/WalkerandNealSimStudyUniformRestrictedRange_Revisedv1AdaptiveSeed81.RData")
# Find the Walker Improvements 
Walkerl1 <- as_tibble(Improvement[,1,])
Walkerl1$loss <- "l1"

Walkerl2 <- as_tibble(Improvement[,2,])
Walkerl2$loss <- "l2"
WalkerComb <- rbind(Walkerl1, Walkerl2)

WalkerLong <- pivot_longer(WalkerComb,
                           cols = -"loss",
                           names_to = "Size",
                           values_to = "Improvement")
WalkerLong$type <- "Walker Slice DPMM"

# Find the Neal Improvements
Improvement <- 100 * (1 - (NPLoss/IndCalLoss))
Neall1 <- as_tibble(Improvement[,1,])
Neall1$loss <- "l1"

Neall2 <- as_tibble(Improvement[,2,])
Neall2$loss <- "l2"
NealComb <- rbind(Neall1, Neall2)

NealLong <- pivot_longer(NealComb,
                         cols = -"loss",
                         names_to = "Size",
                         values_to = "Improvement")
NealLong$type <- "Neal P\u{F2}lya Urn"

BothMethodsLong <- rbind(WalkerLong, NealLong)

loss.labs <- c(TeX("$l_1$-loss"), TeX("$l_2$-loss"))
loss.labs <- c("l[1]-loss", "l[2]-loss")
names(loss.labs) <- c("l1", "l2")

scalefac <- 1.6
pdf("UniformSimStudy.pdf", width = scalefac * 33/4, height = scalefac * 0.25 * 47/4)
punif <- ggplot(BothMethodsLong, aes(x=Size, y=Improvement, fill=type)) +
  geom_boxplot() + facet_wrap(~ loss, labeller = labeller(loss = as_labeller(loss.labs, label_parsed)))
punif <- punif + scale_fill_brewer(palette="RdBu") + scale_x_discrete(limits=c("50", "100", "200", "500"))
punif <- punif + labs(x = TeX("Sample Size $n$"),
                y = "Percentage Improvement \n (Non-Parametric over Independent Calibration)",
                title = "Example 3: Uniform phase",
                subtitle = TeX("$\\theta_i \\in \\left[100, 15000\\right]$"))

punif <- punif + theme(legend.position = "bottom", 
                 plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
punif <- punif + labs(tag = "(c)")
punif <- punif + geom_hline(yintercept = 0, col = "blue")
print(punif)
dev.off()









#########################################################
# Unused but alternative plotting style
l1losses <- BothMethodsLong %>% filter(loss == "l1")

p1 <- ggplot(l1losses, aes(x=Size, y=Improvement, fill=type)) +
  geom_boxplot()
p1 <- p1 + scale_fill_brewer(palette="RdBu") + scale_x_discrete(limits=c("50", "100", "200", "500"))
p1 <- p1 + labs(x = TeX("Sample Size $n$"),
                y = "Improvement over Independent Calibration",
                title = "Three clusters", 
                subtitle = TeX("$l_1$-loss"))
p1 <- p1 + theme(legend.position = "bottom")


l2losses <- BothMethodsLong %>% filter(loss == "l2")
p2 <- ggplot(l2losses, aes(x=Size, y=Improvement, fill=type)) +
  geom_boxplot()
p2 <- p2 + scale_fill_brewer(palette="RdBu") + scale_x_discrete(limits=c("50", "100", "200", "500"))
p2 <- p2 + labs(x = TeX("Sample Size $n$"),
                y = "Improvement over Independent Calibration",
                title = "", 
                subtitle = TeX("$l_2$-loss"))
p2 <- p2 + theme(legend.position = "bottom", legend.title=element_blank(), 
                 legend.text = element_blank(), legend.key = element_rect(fill = "white", colour = "white"))


grid.arrange(p1, p2, nrow=1)


