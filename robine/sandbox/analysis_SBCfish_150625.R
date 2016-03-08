#################################################
# Analyze SBC fish data, test for scale x time interaction
# Author: Robin Elahi
# Date: 150625
#################################################

library(ggplot2)
library(nlme)
library(arm)
library(AICcmodavg)
library(dplyr)
options(dplyr.print_max = 1e9)

rm(list=ls(all=TRUE)) # removes all previous material from R's memory

# get data for analysis
dat <- read.csv("./data/btsDF_sbcfish.csv", 
                header=TRUE, na.strings= c("NA", "NULL"))

dat$dateR <- as.Date(dat$date, format = "%Y-%m-%d")

# make sure the data are ordered by time within time-series (
# necessary for autocorrelation 
head(dat$dateR)

with(dat, table(st, YEAR))


# quick plots to examine data




############################################################
############################################################
###FIGURE 1B: LN(S) V YEAR
############################################################
############################################################
# FIGURE 1B: ln(S) v year; WITH THE OVERALL SLOPE

# Simplest model - what is the overall slope, accounting for random effects?
lmeDat <- richDat

simple <- lme(fixed = richLN ~ I(year0Z + 1990), 
	data = lmeDat, method = "REML", 
	random =  randList1, 
	correlation = corAR1())
summary(simple)

# check normality and homogeneity of variances
par(mfrow = c(1,2))
qqnorm(resid(simple))
qqline(resid(simple))
plot(resid(simple) ~ fitted(simple)); abline(h=0)
plot(simple)

int <- fixed.effects(simple)[1]
slope <- fixed.effects(simple)[2]

ggDat <- lmeDat
ggDat$simple.fitted <- predict(simple)
ggDat$simple.resids <- with(ggDat, simple.fitted - richLN)
plot(simple.fitted ~ year0Z, ggDat)
plot(simple.resids ~ year0Z, ggDat)
year <- ggDat$year0 + 1962

# 
no_legend <- theme(legend.position = "none")

theme1 <- theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 14)) 


no_grid <- theme(panel.grid.major = element_blank()) + 
  theme(panel.grid.minor = element_blank())

# Points colored by studyName, but fitted lines by subSiteID
# With abundance, solid lines; w/o = dashed lines

richLNp2 <- ggplot(ggDat, aes(year0 + 1962, richLN)) + theme1 + no_legend + 
  geom_point(aes(color = studyName), size = 1, alpha = 0.8) + 
   geom_smooth(data = subset(ggDat, AbundYes == 1), 
              method = "lm", se = F, size = 0.15, 
              aes(group = subSiteID, color = studyName)) +
  geom_smooth(data = subset(ggDat, AbundYes == 0), 
              method = "lm", se = F, size = 0.15, linetype = "dashed", 
              aes(group = subSiteID, color = studyName))

ULClabel <- theme(plot.title = element_text(hjust = -0.05, vjust = 1, 
	size = rel(1.5)))

Fig1B <- richLNp2 + xlab("Year") + ylab("ln (S)") +
  labs(title = "B") + ULClabel + 
  geom_abline(intercept = int, slope = slope, color = "black", size = 0.3) + no_grid

# Fig1B
# multiplot(map1, Fig1B, mocFig, layout = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))

#######################################################
# Get model averaged coefficients for 3 datasets, Figure 1C
#########################
# Richness, full dataset
#######################################################

### best model set, CumWt = 0.9

lmeDat <- richDat
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")

# Don't use Gelman's method because SD = 0.5, but I want SD = 1 
# (more similar to the observed SDs for initialRichLN, driver; but not trophic and scale...)

# Rescale duration to visualize coefficient better
# lmeDat$RS.durationZ <- rescale(lmeDat$duration)
lmeDat$RS.durationZ <- scale(lmeDat$duration)[,1]

head(scale(lmeDat$duration))
scale(lmeDat$duration)[,1]

# Rescale year to match duration
# lmeDat$RS.year0Z <- rescale(lmeDat$year0Z)
lmeDat$RS.year0Z <- scale(lmeDat$year0)[,1]

summary(lmeDat)

# How do the standardized inputs compare?
lmeDat %>% summarise(richLNSD = sd(richLN), 
                     initialRichLNSD = sd(initialRichLN), 
                     yearSD = sd(RS.year0Z), 
                     durationSD = sd(RS.durationZ), 
                     driverSD = sd(Prediction),
                     trophicSD = sd(Trophic), 
                     scaleSD = sd(Scale))

# Global model
cm1S <- lme(fixed = richLN ~ 1 +
              RS.year0Z*Scale + RS.year0Z*Prediction +
              RS.year0Z*Trophic + RS.year0Z*initialRichLN + 
              RS.year0Z*RS.durationZ, 
	data = lmeDat, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
summary(cm1S)$tTable
plot(cm1S)

# cm2S <- update(cm1S, ~ . - year0Z:Scale)
cm4S <- update(cm1S, ~ . - RS.year0Z:Trophic)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)

bestModSet <- list(cm1S,cm4S, cm6S)
bestModSet
bestNames <- paste("bestMod", c(1, 4, 6), sep = " ")
bestNames

#######################################################
# write a loop to extract the rowname, coefficient, and upper and lower CL
fullModel <- cm1S

# save the rownames
coefNames <- names(fixed.effects(fullModel))
coefNames  
length(coefNames)
AllReps <- coefNames[11:18] # include only interactions
AllReps #
N <- length(AllReps); N

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)

bestModSet
# get value for year0z
yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, "RS.year0Z", 
                    exclude = list(AllReps), second.ord = FALSE)
yearZcoef
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
yearZname <- "year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative","Year * Driver-Neutral", "Year * Driver-Positive", "Year * Trophic Level 2", "Year * Trophic Level 3", "Year * Initial ln(S)", "Year * Duration")
loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)
loopDat4$ci <- as.numeric(loopDat4$ci)

#######################################################
# save fullData and reduced data
fullBetas <- loopDat4
fullBetas
fullBetas2 <- fullBetas[with(fullBetas, order(modAvgBeta)), ]
fullBetas2
#######################################################
### PREPARE PLOT FOR PPT
#######################################################
# this will reorder based on fullBetas above
loopDat4$xlab2 <- factor(loopDat4$xlabels, 
                         levels = rev(c("Year * Driver-Negative", "Year * Initial ln(S)", 
                                        "Year * Trophic Level 3", "Year * Duration", 
                                        "Year * Trophic Level 2", "Year * Driver-Neutral", 
                                        "Year * Driver-Positive", "Year * Scale-Gamma", 
                                        "Year" )))
loopDat4

# plot
yAxis1 <- theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14))
dodge <- position_dodge(width = -0.7)
plotSpecs <- theme_bw() + yAxis1 + xAxis1 

p1 <- ggplot(loopDat4, aes(x = xlab2, y = modAvgBeta)) +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") + plotSpecs + 
  geom_errorbar(width = 0, aes(ymin = lowerCL, ymax = upperCL), 
                position = dodge, size = 0.8) +
  geom_point(shape = 21, size = 3, fill = "black") +		
  coord_flip() + ylab("Model-Averaged Coefficients")
p1


#######################################################
# Get model averaged coefficients for 3 datasets, Figure 1C
#########################
# Richness, reduced dataset with abundance
#######################################################

lmeDat <- subDat
unique(lmeDat$Prediction)
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")

# Rescale duration to visualize coefficient better
lmeDat$RS.durationZ <- scale(lmeDat$duration)[,1]

# Rescale year to match duration
lmeDat$RS.year0Z <- scale(lmeDat$year0)[,1]


cm1S <- lme(fixed = richLN ~ 1 +  
              RS.year0Z*Scale + RS.year0Z*Prediction +
              RS.year0Z*Trophic + RS.year0Z*initialRichLN + 
              RS.year0Z*RS.durationZ + abundCS,  
	data = lmeDat, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
fixed.effects(cm1S)
summary(cm1S)

cm4S <- update(cm1S, ~ . - RS.year0Z:Trophic)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)


bestModSet <- list(cm1S, cm4S, cm6S)
bestModSet
bestNames <- paste("bestMod", c(1, 4, 6), sep = " ")
bestNames

#######################################################
# write a loop to extract the rowname, coefficient, and upper and lower CL
fullModel <- cm1S

# save the rownames
coefNames <- names(fixed.effects(fullModel))
coefNames  
length(coefNames)
AllReps <- coefNames[12:19] # include only interactions
AllReps #
N <- length(AllReps); N

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)

bestModSet
# get value for year0z
yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, 
                    "RS.year0Z", exclude = list(AllReps), second.ord = FALSE)
yearZcoef
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
yearZname <- "year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative", "Year * Driver-Neutral", 
             "Year * Driver-Positive", "Year * Trophic Level 2", "Year * Trophic Level 3", 
             "Year * Initial ln(S)", "Year * Duration")

loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)

# plot loopDat quickly
yAxis1 <- theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 12))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12))

ggplot(loopDat4[1:9,], aes(x = xlabels, y = modAvgBeta, group = 1)) +
  ylab("Model-Averaged Coefficients") +
  theme_bw() + yAxis1 + xAxis1 +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(width = 0, aes(ymin = modAvgBeta - ci, ymax = modAvgBeta + ci)) +
  geom_point(shape = 21, size = 3, fill = "black") +
  coord_flip() 

subBetasAbund <- loopDat4

#######################################################
# Get model averaged coefficients for 3 datasets, Figure 1C
#########################
# Richness, reduced dataset without abundance
#######################################################

lmeDat <- subDat
unique(lmeDat$Prediction)
# relevel prediction
lmeDat$Prediction <- relevel(lmeDat$Prediction, ref = "none")
dim(lmeDat)

# Rescale duration to visualize coefficient better
lmeDat$RS.durationZ <- scale(lmeDat$duration)[,1]

# Rescale year to match duration
lmeDat$RS.year0Z <- scale(lmeDat$year0)[,1]

cm1S <- lme(fixed = richLN ~ RS.year0Z +  
	RS.year0Z*Scale + RS.year0Z*Prediction +
	RS.year0Z*Trophic + RS.year0Z*initialRichLN + 
	RS.year0Z*RS.durationZ, 
	data = lmeDat, method = "REML", 
	random =  list(rand1, rand2, rand4), 
	correlation = corAR1())
fixed.effects(cm1S)
summary(cm1S)$tTable

# cm2S <- update(cm1S, ~ . - RS.year0Z:Scale)
# cm4S <- update(cm1S, ~ . - RS.year0Z:Trophic)
cm6S <- update(cm1S, ~ . - RS.year0Z:RS.durationZ)
summary(cm6S)$tTable

bestModSet <- list(cm1S, cm6S)
bestModSet
bestNames <- paste("bestMod", c(1, 6), sep = " ")
bestNames

#######################################################
# write a loop to extract the rowname, coefficient, and upper and lower CL
fullModel <- cm1S

# save the rownames
coefNames <- names(fixed.effects(fullModel))
coefNames  
length(coefNames)
AllReps <- coefNames[11:18] # include only interactions
AllReps #
N <- length(AllReps); N

# one matrix for continuous variables; another matrix for categorical
# continuous
mat1 <- matrix(nrow = N, ncol = 4)
colnames(mat1) <- c("modAvgBeta", "lowerCL", "upperCL", "ci")
head(mat1)

bestModSet
# get value for RS.year0Z
yearZcoef <- modavg(cand.set = bestModSet, modnames = bestNames, 
                    "RS.year0Z", exclude = list(AllReps), second.ord = FALSE)
yearZcoef
yearZbeta <- yearZcoef$Mod.avg.beta
yearZlower <- yearZcoef$Lower.CL
yearZupper <- yearZcoef$Upper.CL
yearZci <- abs(yearZlower - yearZupper)/2
str(yearZci)
yearZname <- "RS.year0Z"
yearZlab <- "Year"
str(yearZupper)

yearZbeta - yearZlower
yearZupper - yearZbeta 
yearZci

yearZrow <- as.data.frame(cbind(yearZname, yearZbeta, yearZlower, yearZupper, yearZci, yearZlab))
names(yearZrow) <- c("coefName", "modAvgBeta", "lowerCL", "upperCL", "ci", "xlabels")
yearZrow[, 2:5] <- as.numeric(as.character(unlist(yearZrow[, 2:5])))
str(yearZrow)
yearZrow

# FOR CONTINUOUS DATA
###
for(i in 1:N) {
  coef.i <- AllReps[i]
  coefAvg <- modavg(cand.set = bestModSet, modnames = bestNames, parm = coef.i, 
                    second.ord = FALSE)
  beta <- coefAvg$Mod.avg.beta
  lower <- coefAvg$Lower.CL
  upper <- coefAvg$Upper.CL
  ci <- abs(lower - upper)/2
  # populate matrix with continuous values
  mat1[i,] <- c(beta, lower, upper, ci) 
  loopDat <- as.data.frame(mat1)
  coefName <- AllReps
  loopDat2 <- cbind(coefName, loopDat)
}
loopDat2

##############
xlabels <- c("Year * Scale-Gamma", "Year * Driver-Negative", "Year * Driver-Neutral", 
             "Year * Driver-Positive", "Year * Trophic Level 2", "Year * Trophic Level 3", 
             "Year * Initial ln(S)", "Year * Duration")
loopDat3 <- cbind(loopDat2, xlabels)
loopDat3
loopDat4 <- rbind(yearZrow, loopDat3)
loopDat4
str(loopDat4)

# plot loopDat quickly
yAxis1 <- theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 12))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 12))

ggplot(loopDat4[1:9,], aes(x = xlabels, y = modAvgBeta, group = 1)) +
  ylab("Model-Averaged Coefficients") +
  theme_bw() + yAxis1 + xAxis1 +
  geom_hline(yintercept = 0, color = "red", linetype = "dashed") +
  geom_errorbar(width = 0, aes(ymin = modAvgBeta - ci, ymax = modAvgBeta + ci)) +
  geom_point(shape = 21, size = 3, fill = "black") +
  coord_flip() 

subBetas <- loopDat4

############################################################
############################################################
### FIGURE 1C
############################################################
############################################################

fullBetas
fullBetas[order(fullBetas$modAvgBeta), ]

subBetasAbund
subBetas

# use these values to scale the size of the points
dim(richDat)
dim(subDat)

Betas <- rbind(fullBetas, subBetas, subBetasAbund)
Betas
Betas$Dataset <- c(rep("Full", 9), rep("Reduced", 18))
Betas$Models <- c(rep("No abundance", 18), rep("With abundance", 9))
Betas$DM <- c(rep("Full", 9), rep("Reduced", 9), rep("Reduced, abundance predictor", 9))

# this will reorder based on fullBetas above
Betas$xlab2 <- factor(loopDat4$xlabels, 
                      levels = rev(c("Year * Driver-Negative", "Year * Initial ln(S)", 
                                     "Year * Duration", "Year * Driver-Neutral", 
                                     "Year * Trophic Level 3", "Year * Driver-Positive", 
                                     "Year * Trophic Level 2", 
                                     "Year * Scale-Gamma", 
                                     "Year" )))

Betas <- droplevels(Betas)

#write.csv(Betas, "Betas.csv")

head(Betas)

# plot
theme1 <- theme_bw(base_size = 16) + 
  theme(axis.text = element_text(size = 14)) 

yAxis1 <- theme(axis.title.y = element_blank(), axis.text.y = element_text(size = 14))
xAxis1 <- theme(axis.title.x = element_text(size = 16), axis.text.x = element_text(size = 14))
dodge <- position_dodge(width = -0.7)
plotSpecs <- theme1 + yAxis1 + xAxis1 

allBetas <- ggplot(Betas, aes(x = xlab2, y = modAvgBeta, group = DM, color = Models, size = Dataset)) +
  geom_hline(yintercept = 0, color = "gray", linetype = "dashed") + plotSpecs + 
  geom_point(shape = 19, position = dodge) + 
  scale_size_manual(values = c(5.073*0.6, 3.477*0.6)) +
  geom_errorbar(width = 0, aes(ymin = lowerCL, ymax = upperCL), position = dodge, size = 0.8) + 
  scale_color_manual(values = c("black", "darkgray")) +
  coord_flip() + ylab("Standardized coefficients") +
  no_grid

ULClabel <- theme(plot.title = element_text(hjust = -0.1, vjust = 1, size = rel(1.5)))

mocFig <- allBetas + theme(legend.justification = c(1, 1), legend.position = c(0.52, 0.45), 
                           legend.text = element_text(size = 10), 
                           legend.title = element_text(size = 10)) +
  labs(title = "C") + ULClabel 
  
mocFig

############################################################
############################################################
### FIGURE 1, 3 panels
############################################################
############################################################

multiplot(map1, Fig1B, mocFig, layout = matrix(c(1,1,2,3), nrow = 2, byrow = TRUE))


