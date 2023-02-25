# Code for statistical analysis of data for 'Spread of networked populations is determined by the interplay between dispersal behavior and habitat configuration'
# Authors: Bronwyn Rayfield, Celina B. Baines, Luis J. Gilarranz, and Andrew Gonzalez


####### Required libraries #################################
# Load libraries
require(plyr)
require(boot)
require(lme4)
require(ggplot2)
require(ggpubr)
require(emmeans)
library(simpleboot)


####### Load model and experiment data #####################################

# Set working directory
setwd("...") #the data is in "data for analysis" folder

# Spatial properties of networks and nodes
d_node.properties <- read.table("Node spatial properties.txt", header = T)# node spatial properties
d_network.properties <- read.table("Network spatial properties.txt", header = T)# network spatial properties

# General model simulation data
d_genmodW.network <- read.table("General model_Weighted_Network OccpTime.txt", header = T)
d_genmodB.network <- read.table("General model_Binary_Network OccpTime.txt", header = T)

# F. candida data - model and experiment
d_exp.node <- read.table("Experiment_Node population size.txt", header = T) # experimental data - population size in each node on each observation day
d_mod.node <- read.table("Model_Node population size.txt", header = T)# model data - population size in each node on each day

####### Define functions and miscellaneous values ####

# 1. set colour palette for figures
cbPalette <- c("#63ACBE", "#F5793A", "#A95AA1")

# 2. define functions
## esimates 2.5th (lwr), 50th (fit), and 97.5th (upr) percentiles of a vector
sumBoot <- function(merBoot) {
  return(
    data.frame(fit = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.5, na.rm=TRUE))),
               lwr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.025, na.rm=TRUE))),
               upr = apply(merBoot$t, 2, function(x) as.numeric(quantile(x, probs=.975, na.rm=TRUE)))
    )
  )
}


####### F. candida - Time to node colonization ###########################

# ESTIMATE TIME TO COLONIZATION

## A) MODEL DATA

### 1. Estimate time to colonization for each node
d_mod.timeto1 <- ddply(subset(d_mod.node, N > 0), .(Config, Rep, ModelRep, node), summarize, Clnz.Time = min(Day))

### 2. Add columns to dataset: landscapeID
d_mod.timeto1$topology <- with(d_mod.timeto1, paste(Config, Rep))
d_mod.timeto1$repID <- with(d_mod.timeto1, paste(Config, Rep, ModelRep))

## B) EXPERIMENT DATA

### 1. Estimate time to colonization
d_exp.node.long <- reshape(d_exp.node, 
                           varying = c("Node1", "Node2", "Node3", "Node4", "Node5", "Node6", "Node7", "Node8", "Node9", "Node10"),
                           v.names = "N",
                           timevar = "node",
                           times = c(1:10),
                           direction = "long",
                           drop = c("Date")
)
                           
d_exp.timeto1 <- ddply(subset(d_exp.node.long, N > 0), .(Config, Rep, node), summarize, Clnz.Time = min(Day))

### 2. Add columns to dataset: landscapeID
d_exp.timeto1$landscapeID <- with(d_exp.timeto1, paste(Config, Rep))


# ESTIMATE EFFECT OF SPATIAL PROPERTIES OF NODES ON TIME TO COLONIZATION

## A) MODEL DATA

### 1. Add spatial data to data frame providing time to colonization
d_mod.timeto1 <- join(d_mod.timeto1, d_node.properties, type = "right")

### 2. Plot time to colonization as a function of distance from source and configuration

#### a. Build liner mixed model colonization time as the response, distance from source and configuration as fixed predictors, and topology and repID as random effects
m_mod.clnz.dist <- lmer(log10(Clnz.Time) ~ poly(distToSource,2) * as.factor(Config) + (1|topology) + (1|repID), data = subset(d_mod.timeto1, node != 5))

hist(summary(m_mod.clnz.dist)$residuals)
plot(summary(m_mod.clnz.dist)$residuals ~ fitted(m_mod.clnz.dist))


#### b. Extract predictions +/- 95% CI of colonization time as a function of distance from source
nd_mod.clnz.dist <- data.frame(unique(subset(d_mod.timeto1, node != 5, select = c(distToSource, Config))), "topology" = NA)

mySumm_mod.clnz.dist <- function(.) {
  predict(., newdata=nd_mod.clnz.dist, re.form=~0)
}

boot_mod.clnz.dist <- bootMer(m_mod.clnz.dist, mySumm_mod.clnz.dist, nsim = 1000, use.u=F, type="parametric")
est_mod.clnz.dist <- cbind(nd_mod.clnz.dist, sumBoot(boot_mod.clnz.dist))

#### c. Untransform fit and confidence interval estimates to response units (days)
est_mod.clnz.dist$fit.untr <- 10^est_mod.clnz.dist$fit
est_mod.clnz.dist$lwr.untr <- 10^est_mod.clnz.dist$lwr
est_mod.clnz.dist$upr.untr <- 10^est_mod.clnz.dist$upr

### 3. Plot colonization time as a function of distance from source
g_mod.clnz.dist <- ggplot()+
  geom_jitter(data = subset(d_mod.timeto1, node != 5), aes(x = distToSource, y = Clnz.Time, col = as.factor(Config)), size = 2, pch = 1)+
  geom_line(data = est_mod.clnz.dist, aes(x = distToSource, y = fit.untr, col = as.factor(Config), lty = as.factor(Config)), size = 2)+
  geom_ribbon(data = est_mod.clnz.dist, aes(x = distToSource, ymin = lwr.untr, ymax = upr.untr, fill = as.factor(Config)), alpha = 0.5)+
  scale_y_continuous("Days to node\ncolonization", limits = c(0, 200), breaks = c(0, 50, 100, 150, 200))+
  scale_x_continuous("Distance from source (mm)", limits = c(0, 500), breaks = c(0, 125, 250, 375, 500))+
  scale_color_manual(guide = F, "Habitat configuration", values = cbPalette, labels = c("Lattice", "Partially random", "Random"))+
  scale_fill_manual(guide = F, "Habitat configuration", values = cbPalette, labels = c("Lattice", "Partially random", "Random"))+
  scale_linetype_manual(guide = F, "Habitat configuration", labels = c("Lattice", "Partially random", "Random"), values = c(2,1,3))+
  labs(tag = "A")+
  theme_bw(base_size=12)+ 
  theme(legend.key.width = unit(1, "cm"), plot.caption = element_text(hjust = 0), strip.background = element_rect(fill = "white", colour = "black"), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))


## B) EXPERIMENT DATA

### 1. Add spatial data to data frame providing time to colonization
d_exp.timeto1 <- join(d_exp.timeto1, d_node.properties)

### 2. Build linear mixed model with colonization time as the response, and distance from source, configuration, and their two-way interaction, and landscapeID as a random effect
m_exp.clnz.dist.intx <- lmer(Clnz.Time ~ distToSource * as.factor(Config) + (1|landscapeID), data = subset(d_exp.timeto1, node != 5))

hist(summary(m_exp.clnz.dist.intx)$residuals)
plot(summary(m_exp.clnz.dist.intx)$residuals ~ fitted(m_exp.clnz.dist.intx))


# linear regression with configuration and distance from source as predictors
m_exp.clnz.dist.add <- lmer(Clnz.Time ~ distToSource + as.factor(Config) + (1|landscapeID), data = subset(d_exp.timeto1, node != 5))

hist(summary(m_exp.clnz.dist.add)$residuals)
plot(summary(m_exp.clnz.dist.add)$residuals ~ fitted(m_exp.clnz.dist.add))


### 3. Extract results of linear mixed models
result_exp.clnz.dist.intx <- drop1(m_exp.clnz.dist.intx, test = "Chisq")
result_exp.clnz.dist.add <- drop1(m_exp.clnz.dist.add, test = "Chisq")


#### 4. Extract predicted mean +/- 95% CI of colonization time as a function of distance from source
nd_exp.clnz.dist <- data.frame(unique(subset(d_exp.timeto1, node != 5, select = c(distToSource, Config))), "landscapeID" = NA)

mySumm_exp.clnz.dist <- function(.) {
  predict(., newdata=nd_exp.clnz.dist, re.form=~0)
}

boot_exp.clnz.dist <- bootMer(m_exp.clnz.dist.add, mySumm_exp.clnz.dist, nsim = 1000, use.u=F, type="parametric")

est_exp.clnz.dist <- cbind(nd_exp.clnz.dist, sumBoot(boot_exp.clnz.dist))

### 5. Plot colonization time as a function of distance from source

g_exp.clnz.dist <- ggplot()+
  geom_ribbon(data = est_exp.clnz.dist, aes(x = distToSource, ymin = lwr, ymax = upr, fill = as.factor(Config)), alpha = 0.3)+
  geom_jitter(data = subset(d_exp.timeto1, distToSource != 0), aes(x = distToSource, y = Clnz.Time, col = as.factor(Config)), size = 2, pch = 1)+
  geom_line(data = est_exp.clnz.dist, aes(x = distToSource, y = fit, col = as.factor(Config), lty = as.factor(Config)), size = 2)+
  scale_y_continuous("Days to node\ncolonization", limits = c(0, 80), breaks = c(0, 20, 40, 60, 80))+
  scale_x_continuous("Distance from source (mm)", limits = c(0, 500), breaks = c(0, 125, 250, 375, 500))+
  scale_color_manual("Habitat configuration", values = cbPalette, labels = c("Lattice", "Partially random", "Random"))+
  scale_fill_manual("Habitat configuration", values = cbPalette, labels = c("Lattice", "Partially random", "Random"))+
  scale_linetype_manual("Habitat configuration", labels = c("Lattice", "Partially random", "Random"), values = c(2,1,3))+
  labs(tag = "B")+
  theme_bw(base_size=12)+ 
  theme(legend.position = c(0.78, 0.85), legend.title = element_text(size = 8), legend.text = element_text(size = 8), legend.key.width = unit(1, "cm"), legend.key.height = unit(0.3, "cm"), plot.caption = element_text(hjust = 0), strip.background = element_rect(fill = "white", colour = "black"), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))


####### F. candida - Time to full network occupancy ##################

# ESTIMATE TIME TO FULL OCCUPANCY FOR EACH NETWORK

## A) MODEL DATA

### 1. Estimate time to full occupancy for each network

#### a. Calculate proportion of nodes that are occupied in each network
d_mod.node$occupied <- with(d_mod.node, ifelse(N > 0, 1, 0))
d_mod.occupancy <- ddply(d_mod.node, .(Config, Rep, ModelRep, Day), summarize, prop.occupancy = mean(occupied))

#### b. Estimate minimum day on which network occupancy = 1
d_mod.spreadtime <- ddply(subset(d_mod.occupancy, prop.occupancy == 1), .(Config, Rep, ModelRep), summarize, Occp.Time = min(Day))
d_mod.spreadtime$topology <- with(d_mod.spreadtime, paste(Config, Rep))

## B) EXPERIMENT DATA

### 1. Estimate time to full occupancy for each network

d_exp.node$full.occp <- with(d_exp.node, ifelse(Node1 > 0 & Node2 > 0 & Node3 > 0 & Node4 >0 & Node5 > 0 & Node6 >0 & Node7 > 0 & Node8 > 0 & Node9 > 0 & Node10 > 0, 1, 0)) #if network is fully occupied, full.occp = 1, 0 otherwise

#### b. Estimate minimum day on which network occupancy = 1

d_exp.spreadtime <- ddply(subset(d_exp.node, full.occp > 0), .(Config, Rep), summarize, Occp.Time = min(Day))
d_exp.spreadtime$landscapeID <- with(d_exp.spreadtime, paste(Config, Rep))


# ESTIMATE EFFECT OF LANDSCAPE CONFIGURATION ON TIME TO FULL OCCUPANCY

## A) MODEL DATA

## SMW and RND configurations

### 1. Build a linear mixed model with time to full occupancy as the response, landscape configuration as a fixed predictor, and topology as a random effect; exclude lattice networks because they have only one topology
m_mod.occp.LC.23 <- lmer(formula = log10(Occp.Time) ~ as.factor(Config) + (1|topology), data = subset(d_mod.spreadtime, Config != 1))

hist(summary(m_mod.occp.LC.23)$residuals)
plot(summary(m_mod.occp.LC.23)$residuals ~ fitted(m_mod.occp.LC.23))


### 2. Extract mean +/- 95% CI time to full occupancy for each landscape configuration (SMW and RND only)
nd_mod.occp.LC.23 <- data.frame("Config" = c(2,3), "topology" = NA)

mySumm_mod.occup.LC <- function(.) {
  predict(., newdata=nd_mod.occp.LC.23, re.form=~0)
}

est_mod.occp.LC.23 <- cbind(nd_mod.occp.LC.23, sumBoot(bootMer(m_mod.occp.LC.23, mySumm_mod.occup.LC, nsim = 1000, use.u=F, type="parametric")))


## LAT configuration

### 1. Build a linear intercept model with time to full occupancy as the response
d_mod.occp.LC.1 <- subset(d_mod.spreadtime, Config == 1) # subset dataframe to include lattice networks only
d_mod.occp.LC.1$log10_Occp.Time <- with(d_mod.occp.LC.1, log10(Occp.Time)) # create column - log10(Occp.Time)

m_mod.occp.LC.1_simple <- lm(formula = log10_Occp.Time ~ 1, data = d_mod.occp.LC.1) # lm.boot cannot handle functions (like log10) within the formula; this is the same as m_mod.occp.LC.1 above, but with functions on variables performed outside of/before the model

mboot.mod.occp.LC.1 <- lm.boot(m_mod.occp.LC.1_simple, rows = F, R = 1000)
mboot.fits_mod.occp.LC.1 <- data.frame("coef" = samples(mboot.mod.occp.LC.1, name = c("coef")))

est_mod.occp.LC.1 <- data.frame("Config" = 1, "fit" = mean(mboot.fits_mod.occp.LC.1$coef), "lwr" = quantile(mboot.fits_mod.occp.LC.1$coef, 0.025), "upr" = quantile(mboot.fits_mod.occp.LC.1$coef, 0.975))


## Combine mean +/- 95% confidence level estimates of time to full occupancy for each landscape configuration into one dataset
est_mod.occp.LC <- rbind.fill(est_mod.occp.LC.1, est_mod.occp.LC.23)

### 3. Untransform the mean +/- 95% CI values to original response units
est_mod.occp.LC$fit.untr <- 10^est_mod.occp.LC$fit
est_mod.occp.LC$lwr.untr <- 10^est_mod.occp.LC$lwr
est_mod.occp.LC$upr.untr <- 10^est_mod.occp.LC$upr



### 4. Plot time to full occupancy as a function of landscape configuration
g_mod.occp.LC <- ggplot(est_mod.occp.LC, aes(x = as.factor(Config), y = fit.untr, fill = as.factor(Config)))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lwr.untr, ymax = upr.untr), width = 0.1, size = 0.5)+
  scale_y_continuous("Days to full network\noccupancy", limits = c(0, 80), breaks = c(0, 20, 40, 60, 80))+
  scale_x_discrete("", labels = c("Lattice", "Partially random", "Random"))+
  scale_fill_manual(guide = F, values = cbPalette)+
  labs(tag = "A")+
  theme_bw(base_size=12)+ 
  theme(plot.caption = element_text(hjust = 0), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))

### 5. Estimate effect size of landscape configuration on time to full occupancy
ef_mod.occp.LC_LATSMW <- est_mod.occp.LC[2,8] / est_mod.occp.LC[1,8] # Lattice vs Partially random
ef_mod.occp.LC_LATRND <- est_mod.occp.LC[3,8] / est_mod.occp.LC[1,8] # Lattice vs Random
ef_mod.occp.LC_SMWRND <- est_mod.occp.LC[3,8] / est_mod.occp.LC[2,8] # Partially random vs Random

## B) EXPERIMENT DATA

### 1. Estimate effect of landscape configuration on time to full occupancy
m_exp.occp.LC <- lm(formula = log10(Occp.Time) ~ as.factor(Config), data = d_exp.spreadtime)
summary_exp.occp.LC <- summary(m_exp.occp.LC)


### 2. Extract results of linear model
result_exp.occp.LC <- anova(m_exp.occp.LC)

### 3. Estimate mean time to full occupancy +/- 95% confidence intervals for each landscape configuration
est_exp.occp.LC <- as.data.frame(emmeans(m_exp.occp.LC, "Config"), type = "response")

### 4. Plot time to full occupancy as a function of landscape configuration

g_exp.occp.LC <- ggplot(est_exp.occp.LC, aes(x = as.factor(Config), y = response, fill = as.factor(Config)))+
  geom_bar(stat = "identity")+
  geom_errorbar(aes(ymin = lower.CL, ymax = upper.CL), width = 0.1, size = 0.5)+
  #geom_image(data = folsomia_image.exp.occp.LC, aes(image = image, x = Config, y = y), size = 0.2) +
  #facet_wrap(~type)+
  scale_y_continuous("", limits = c(0, 80), breaks = c(0, 20, 40, 60, 80))+
  scale_x_discrete("", labels = c("Lattice", "Partially random", "Random"))+
  scale_fill_manual(guide = F, values = cbPalette)+
  labs(tag = "B")+
  theme_bw(base_size=12)+ 
  theme(plot.caption = element_text(hjust = 0), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))

### 5. Estimate contrasts of each configuration against each other configuration
con_exp.occp.LC <- contrast(regrid(emmeans(m_exp.occp.LC, "Config")), 'tukey')

### 6. Estimate effect size of configuration on time to full occupancy
ef_exp.occp.LC_LATSMW <- est_exp.occp.LC[2,2] / est_exp.occp.LC[1,2] # Lattice vs Partially random
ef_exp.occp.LC_LATRND <- est_exp.occp.LC[3,2] / est_exp.occp.LC[1,2] # Lattice vs Random
ef_exp.occp.LC_SMWRND <- est_exp.occp.LC[3,2] / est_exp.occp.LC[2,2] # Partially random vs Random


# ESTIMATE EFFECT OF NETWORK SPATIAL PROPERTIES ON TIME TO FULL OCCUPANCY

### Join time to full occupancy data with network properties data
d_mod.spreadtime <- join(d_mod.spreadtime, d_network.properties)
d_exp.spreadtime <- join(d_exp.spreadtime, d_network.properties)

# Algebraic connectivity (weighted)
## A) MODEL DATA


### 1. Build linear mixed model with occupancy time as the response, algebraic connectivity as a fixed predictor, and topology as a random effect
m_mod.occp.AC <- lmer(log10(Occp.Time) ~ log10(AlgebraicConnectivity_Weighted) + (1|topology), data = d_mod.spreadtime)

hist(summary(m_mod.occp.AC)$residuals)
plot(summary(m_mod.occp.AC)$residuals ~ fitted(m_mod.occp.AC))


### 2. Extract predicted mean +/- 95% CI of occupancy time as a function of algebraic connectivity
nd_mod.occp.AC <- data.frame("AlgebraicConnectivity_Weighted" = unique(d_network.properties$AlgebraicConnectivity_Weighted), "topology" = NA)

mySumm_mod.occp.AC <- function(.) {
  predict(., newdata=nd_mod.occp.AC, re.form=~0)
}

boot_mod.occp.AC <- bootMer(m_mod.occp.AC, mySumm_mod.occp.AC, nsim = 1000, use.u=F, type="parametric")
est_mod.occp.AC <- cbind(nd_mod.occp.AC, sumBoot(boot_mod.occp.AC))

### 3. Untransform predicted values to response units (days)
est_mod.occp.AC$fit_untr <- 10^est_mod.occp.AC$fit
est_mod.occp.AC$lwr_untr <- 10^est_mod.occp.AC$lwr
est_mod.occp.AC$upr_untr <- 10^est_mod.occp.AC$upr

### 4. Plot occupancy time as a function of algebraic connectivity
g_mod.occp.AC <- ggplot()+
  geom_point(data = d_mod.spreadtime, aes(x = AlgebraicConnectivity_Weighted, y = Occp.Time, col = as.factor(Config)), size = 2, pch = 1)+
  geom_ribbon(data = est_mod.occp.AC, aes(x = AlgebraicConnectivity_Weighted, ymin = lwr_untr, ymax = upr_untr), alpha = 0.3)+
  geom_line(data = est_mod.occp.AC, aes(x = AlgebraicConnectivity_Weighted, y = fit_untr), size = 1.5)+
  scale_y_log10("Days to full network\noccupancy", limits = c(8,200), breaks = c(10, 50, 100, 200))+
  scale_x_continuous("Algebraic connectivity")+
  scale_color_manual("Habitat configuration", labels = c("Lattice", "Partially random", "Random"), values = cbPalette)+
  labs(tag = "C")+
  theme_bw(base_size=12)+ 
  theme(legend.position = c(0.7, 0.8), plot.caption = element_text(hjust = 0), strip.text = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))


## B) EXPERIMENT DATA

### 1. Build linear model with occupancy time as the response and algebraic connectivity as a fixed predictor
m_exp.occp.AC <- lm(log10(Occp.Time) ~ AlgebraicConnectivity_Weighted, data = d_exp.spreadtime)

hist(summary(m_exp.occp.AC)$residuals)
plot(summary(m_exp.occp.AC)$residuals ~ fitted(m_exp.occp.AC))

### 2. Extract results of linear model
result_exp.occp.AC <- anova(m_exp.occp.AC)
summary_exp.occp.AC <- summary(m_exp.occp.AC)

### 3. Extract predicted mean +/- 95% CI of occupancy time as a function of algebraic connectivity
d_exp.spreadtime$log10_Occp.Time <- with(d_exp.spreadtime, log10(Occp.Time))
m_exp.occp.AC_simple <- lm(log10_Occp.Time ~ AlgebraicConnectivity_Weighted, data = d_exp.spreadtime) # lm.boot can't handle functions within the formulas; rewriting this model to simplify it

xpts_AC <- data.frame("AlgebraicConnectivity_Weighted" = unique(d_network.properties$AlgebraicConnectivity_Weighted))# create dataframe of algebraic connectivity values at which response will be predicted
mboot_exp.occp.AC <- lm.boot(m_exp.occp.AC_simple, R = 1000, rows = T, new.xpts = xpts_AC) # bootstrap the linear model by resampling rows of the original data
mboot.fits_exp.occp.AC <- fitted(mboot_exp.occp.AC) # extract bootstrap replicate values for each value of AC

## create dataframe with estimates of fit (mean of the boostrap replicates for each value of AC), lower 95% confidence limit (2.5th percentile of bootstrap replicates) and upper 95% confidence limit (97.5th percentile of boostrap replicates)
est_exp.occp.AC <- NULL
for(i in 1:nrow(xpts_AC)){
  fit <- mean(mboot.fits_exp.occp.AC[i,])
  lwr <- quantile(mboot.fits_exp.occp.AC[i,], 0.025)
  upr <- quantile(mboot.fits_exp.occp.AC[i,], 0.975)
  est_exp.occp.AC <- rbind(est_exp.occp.AC, data.frame("AlgebraicConnectivity_Weighted" = xpts_AC$AlgebraicConnectivity_Weighted[i], "fit" = fit, "lwr" = lwr, "upr" = upr))
}

### 4. Untransform predicted values to response units (days)

est_exp.occp.AC$fit.untr <- 10^est_exp.occp.AC$fit
est_exp.occp.AC$lwr.untr <- 10^est_exp.occp.AC$lwr
est_exp.occp.AC$upr.untr <- 10^est_exp.occp.AC$upr



### 5. Plot occupancy time as a function of algebraic connectivity

g_exp.occp.AC <- ggplot()+
  geom_point(data = d_exp.spreadtime, aes(y = Occp.Time, x = AlgebraicConnectivity_Weighted, col = as.factor(Config)), size = 2, pch = 1)+
  geom_ribbon(data = est_exp.occp.AC, aes(x = AlgebraicConnectivity_Weighted, ymin = lwr.untr, ymax = upr.untr), alpha = 0.3)+
  geom_line(data = est_exp.occp.AC, aes(x = AlgebraicConnectivity_Weighted, y = fit.untr), size = 1.5)+
  scale_y_log10("", limits = c(8,200), breaks = c(10, 50, 100, 200))+
  scale_x_continuous("Algebraic connectivity")+
  scale_color_manual(guide = F, labels = c("Lattice", "Partially random", "Random"), values = cbPalette)+
  labs(tag = "D")+
  theme_bw(base_size=12)+ 
  theme(plot.caption = element_text(hjust = 0), strip.text = element_blank(), strip.background = element_blank(), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))



# Compare spatial metrics

## Experiment
### Binary algebraic connectivity
### 1. Build linear model with occupancy time as the response and algebraic connectivity as a fixed predictor
m_exp.occp.ACB <- lm(log10(Occp.Time) ~ AlgebraicConnectivity_Binary, data = d_exp.spreadtime)

hist(summary(m_exp.occp.ACB)$residuals)
plot(summary(m_exp.occp.ACB)$residuals ~ fitted(m_exp.occp.ACB))

### 2. Extract results of linear model
result_exp.occp.ACB <- anova(m_exp.occp.ACB)
summary_exp.occp.ACB <- summary(m_exp.occp.ACB)

### Binary network diameter
### 1. Build linear model with occupancy time as the response and network diameter as a fixed predictor
m_exp.occp.NDB <- lm(log10(Occp.Time) ~ as.factor(NetworkDiameter_Binary), data = d_exp.spreadtime)

hist(summary(m_exp.occp.NDB)$residuals)
plot(summary(m_exp.occp.NDB)$residuals ~ fitted(m_exp.occp.NDB))

### 2. Extract results of linear model
result_exp.occp.NDB <- anova(m_exp.occp.NDB)
summary_exp.occp.NDB <- summary(m_exp.occp.NDB)

### Weighted network diameter
### 1. Build linear model with occupancy time as the response and network diameter as a fixed predictor
m_exp.occp.ND <- lm(log10(Occp.Time) ~ NetworkDiameter_Weighted, data = d_exp.spreadtime)

hist(summary(m_exp.occp.ND)$residuals)
plot(summary(m_exp.occp.ND)$residuals ~ fitted(m_exp.occp.ND))

### 2. Extract results of linear model
result_exp.occp.ND <- anova(m_exp.occp.ND)
summary_exp.occp.ND <- summary(m_exp.occp.ND)

## Model

mod.spearmanR.ACB <- with(d_mod.spreadtime, cor(Occp.Time, AlgebraicConnectivity_Binary, method = "spearman"))
mod.spearmanR.AC <- with(d_mod.spreadtime, cor(Occp.Time, AlgebraicConnectivity_Weighted, method = "spearman"))
mod.spearmanR.NDB <- with(d_mod.spreadtime, cor(Occp.Time, NetworkDiameter_Binary, method = "spearman"))
mod.spearmanR.ND <- with(d_mod.spreadtime, cor(Occp.Time, NetworkDiameter_Weighted, method = "spearman"))


####### General model analysis - spread (network occupancy time) as a function of B, the exponent of the dispersal kernel #############

# We are removing networks from the stochastic model output that do not reach full occupancy within the timeframe of the simulations (300 days).
# In stochastic model, when full network occupancy does not occur in 300 time steps, Occp.Time is set to 0.
# Convert Occp.Time == 0 to Occp.Time is NA
d_genmodW.network$Occp.Time[d_genmodW.network$Occp.Time == 0] <- NA
d_genmodB.network$Occp.Time[d_genmodB.network$Occp.Time == 0] <- NA


# Estimate mean occupancy time as a function of configuration and B
### 1. Set seed
set.seed(2020)

### 2. Create function - estimate mean of a vector of node population size
fc_Occp.mean <- function(d, i){
  d2 <- d[i,]
  return(mean(d2$Occp.Time, na.rm = T))
}

# WEIGHTED MODEL

# nb: using bootstrapped estimates of time to full occupancy rather than a linear model approach because data do not meet assumptions of linear regression (non-normal residuals even after transformations attempted)

### 1. Use for loop to get bootstrapped estimates of the mean occupancy time for each configuration and value of B
d_genmodW.OccpTimeB <- NULL


for(c in 1:3){
  for(i in unique(d_genmodW.network$B)){
    dSUB <- subset(d_genmodW.network, Config == c & B == i)
    Occp.mean <- boot::boot(dSUB, fc_Occp.mean, R = 1000)$t0
    d_genmodW.OccpTimeB <- rbind(d_genmodW.OccpTimeB, data.frame("Config" = c, "B" = i, "meanOccp" = c(Occp.mean)))
  }
}


## 2. Plot time to full network occupancy as a function of habitat configuration and the exponent of the dispersal kernel

g_genmodW.occpB <- ggplot()+
  # collembola parameter range
  geom_rect(aes(xmin=-0.0430 - 0.007, xmax=-0.0430 + 0.007, ymin=0, ymax=Inf), alpha = 0.2)+ # ymin = 0 introduces infinite values in y-axis
  geom_vline(xintercept = -0.0430, size = 3)+
  # model data
  geom_jitter(data = subset(d_genmodW.network, Config != 1 & B > -0.043 - 0.008), aes(x = B, y = Occp.Time, col = as.factor(Config)), width = 0.0005, height = 0, size = 1, alpha = 0.5)+
  geom_jitter(data = subset(d_genmodW.network, Config == 1 & B > -0.043 - 0.008), aes(x = B, y = Occp.Time, col = as.factor(Config)), width = 0.0005, height = 0, size = 1, alpha = 0.5)+
  geom_line(data = subset(d_genmodW.OccpTimeB, B > -0.043 - 0.008), aes(x = B, y = meanOccp, col = as.factor(Config)), size = 2)+
  # aesthetics
  labs(tag = "B")+
  scale_y_log10("", limits = c(3,300), breaks = c(3, 10, 30, 100, 300))+
  scale_x_continuous(expression(paste("Exponent of the dispersal kernel, ", italic(b))))+
  scale_color_manual("Habitat configuration", labels = c("Lattice", "Partially random", "Random"), values = cbPalette)+
  theme_bw(base_size=12)+ 
  theme(legend.position = c(0.8, 0.8), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))

# BINARY MODEL

### 1. Use for loop to get bootstrapped estimates of the mean occupancy time for each configuration and value of B
d_genmodB.OccpTimeB <- NULL


for(c in 1:3){
  for(i in unique(d_genmodB.network$B)){
    dSUB <- subset(d_genmodB.network, Config == c & B == i)
    Occp.mean <- boot::boot(dSUB, fc_Occp.mean, R = 1000)$t0
    d_genmodB.OccpTimeB <- rbind(d_genmodB.OccpTimeB, data.frame("Config" = c, "B" = i, "meanOccp" = c(Occp.mean)))
  }
}


## 2. Plot time to full network occupancy as a function of habitat configuration and the exponent of the dispersal kernel

g_genmodB.occpB <- ggplot()+
  # collembola parameter range
  geom_rect(aes(xmin=-0.0430 - 0.007, xmax=-0.0430 + 0.007, ymin=0, ymax=Inf), alpha = 0.2)+
  geom_vline(xintercept = -0.0430, size = 3)+
  # model data
  geom_jitter(data = subset(d_genmodB.network, B > -0.043 - 0.008), aes(x = B, y = Occp.Time, col = as.factor(Config)), size = 1, alpha = 0.5, width = 0.0005, height = 0)+
  geom_line(data = subset(d_genmodB.OccpTimeB, B > -0.043 - 0.008), aes(x = B, y = meanOccp, col = as.factor(Config)), size = 2)+
  # aesthetics
  labs(tag = "A")+
  scale_y_log10("Days to full network occupancy", limits = c(3,300), breaks = c(3, 10, 30, 100, 300))+
  scale_x_continuous(expression(paste("Exponent of the dispersal kernel, ", italic(b))))+
  scale_color_manual("Habitat configuration", labels = c("Lattice", "Partially random", "Random"), values = cbPalette)+
  theme_bw(base_size=12)+ 
  theme(legend.position = c(0.8, 0.8), panel.border = element_rect(colour="black"), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"), axis.text = element_text(colour = "black"))



####### Manuscript Results ############

# Main manuscript

# F. candida - effect of landscape configuration on time to full network occupancy
## Model
### Effect size of landscape configuration on time to full network occupancy
ef_mod.occp.LC_LATSMW # Lattice vs Partially random
ef_mod.occp.LC_LATRND # Lattice vs Random
ef_mod.occp.LC_SMWRND # Partially random vs Random

# F. candida - effect of landscape configuration on time to full network occupancy
## Experiment
### LM results of effect of habitat configuration on time to full network occupancy
summary_exp.occp.LC
result_exp.occp.LC
### Effect size of habitat configuration on time to node colonization
ef_exp.occp.LC_LATSMW # Lattice vs Partially random
ef_exp.occp.LC_LATRND # Lattice vs Random
ef_exp.occp.LC_SMWRND # Partially random vs Random
### Pariwise tests - differences between habitat configurations in time to full network occupancy
con_exp.occp.LC

# F. candida - effect of spatial properties on time to full network occupancy
## Experiment
### Effect of weighted algebraic connectivity on time to full network occupancy
summary_exp.occp.AC
result_exp.occp.AC

# Appendix

# F. candida - comparison of effect of spatial properties on time to full network occupancy
## Model
### Effect of binary algebraic connectivity on time to full network occupancy
mod.spearmanR.ACB
### Effect of weighted algebraic connectivity on time to full network occupancy
mod.spearmanR.AC
### Effect of binary network diameter on time to full network occupancy
mod.spearmanR.NDB
### Effect of weighted network diameter on time to full network occupancy
mod.spearmanR.ND

# F. candida - comparison of effect of spatial properties on time to full network occupancy
## Experiment
### Effect of binary algebraic connectivity on time to full network occupancy
result_exp.occp.ACB
summary_exp.occp.ACB
### Effect of weighted algebraic connectivity on time to full network occupancy
result_exp.occp.AC
summary_exp.occp.AC
### Effect of binary network diameter on time to full network occupancy
result_exp.occp.NDB
summary_exp.occp.NDB
### Effect of weighted network diameter on time to full network occupancy
result_exp.occp.ND
summary_exp.occp.ND

# F. candida - effects of spatial properties on time to node colonization
## Experiment
### LMM results of effect of habitat configuration and distance from source on time to node colonization
result_exp.clnz.dist.add


####### Figures ###########

# 1. Set working directory
setwd("...") # folder that you want to save the figures in

# 2. Build figures

## Figure 2

f2 <- ggarrange(g_mod.occp.LC, g_exp.occp.LC, g_mod.occp.AC, g_exp.occp.AC, labels = c(" Model       ", "Experiment", "", ""), label.x = 0.2, font.label = list(size = 10, face = "plain"), common.legend = F, align = "hv")

tiff("Figure 2.tiff", width = 20, height = 20, units = "cm", compression = "none", res = 600)
f2
dev.off()

## Figure 3

f3 <- ggarrange(g_genmodB.occpB, g_genmodW.occpB, labels = c("Binary       ", "Weighted", "", ""), label.x = 0.2, vjust = -0.01, font.label = list(size = 10, face = "plain"), common.legend = T, ncol = 2, align = "hv")

tiff("Figure 3.tiff", width = 22, height = 14, units = "cm", compression = "none", res = 600)
f3
dev.off()

# Appendix

# Days to node colonization as a function of distance from source
fS10 <- ggarrange(g_mod.clnz.dist, g_exp.clnz.dist, labels = c(" Model       ", "Experiment", "", ""), label.x = 0.2, font.label = list(size = 10, face = "plain"), common.legend = F, align = "hv")

tiff("Figure_distance from source.tiff", width = 22, height = 10, units = "cm", compression = "none", res = 600)
fS10
dev.off()

