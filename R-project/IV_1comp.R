#IV_1comp.R
#This script fits PK data following IV administration with a 1-compartment model.
#Written by: David Dubins
#Date: 04-Feb-26

#Read in the data (Time, Concentration) for a Subject
#For PC:
PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

#Replace 300g with mouse weight, or replace the expression below with dose in mg 
Dose <- 15*300/1000 #dose in mg/kg * mouse weight in g/1000 = Dose in mg

# Define the model you are using to fit the data:
iv_1comp <- function(Vd,L1) {
  (Dose/Vd)*exp(-L1*PKexample$Time)  #Dose in mg, Vd in L
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((PKexample$Concentration-iv_1comp(p[1],p[2]))^2)
#We are using "p" as an array of parameters
#p[1] = Vd (in mL)
#p[2] = L1 (in 1/h)

#Plot the original data on a graph:
plot(PKexample$Time,PKexample$Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (µg/mL)")

# Try your hand at fitting your parameters to get in the right ballpark. 
# Myfity has your first guesses at g1, g2, g3, g4, g5, g6, g7 (replace 1, 2, 3, 4, 5, 6, 7)
# with your first guesses)
g1 <- 0.6        #first guess for Vd in L
g2 <- 0.65125    #first guess for L1

# Plot your guess on the plot:
Myfitx <- PKexample$Time
Myfity <- iv_1comp(g1,g2)
lines(spline(Myfitx, Myfity))
#if your first guesses are off the plot, find them:
#plot(Myfitx, Myfity)

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2), hessian = TRUE)
# See the results of your fit:
Myfitout

#Now superimpose our minimized fit on a brand new plot:
plot(PKexample$Time,PKexample$Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (µg/mL)")
PKfitx <- PKexample$Time
PKfity <- iv_1comp(Myfitout$estimate[1],Myfitout$estimate[2])
lines(spline(PKfitx, PKfity))

#Now superimpose our minimized fit on a semi-log plot:
plot(PKexample$Time,log(PKexample$Concentration),main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (µg/mL)")
PKfitx <- PKexample$Time
PKfity_ln <- log(iv_1comp(Myfitout$estimate[1],Myfitout$estimate[2]))
lines(spline(PKfitx, PKfity_ln))

#Calculate the Covariance matrix and the correlation matrix
covmat <- 2 * Myfitout$minimum / (length(PKexample$Concentration) - 2) * solve(Myfitout$hessian)
#Convert Covariance Matrix to the Correlation Matrix.
cormat <- cov2cor(covmat)
param.names <- c("Vd", "L1")
dimnames(cormat) <- list(param.names, param.names)

# Calculate the standard errors of parameter estimates
semat <- sqrt(diag(covmat))
names(semat) <- c("Vd (L)","L1 (1/h)")

# Report the parameter estimates and their standard errors:
results <- data.frame(
  Estimate  = round(as.numeric(Myfitout$estimate),4), #Adjust rounding as needed
  StdError  = round(semat,4)   #second argument is # decimals
)
colnames(results) <- c("Estimate","±StdErr" )
results

#Calculate the rank of the fit, p_eff
p_eff <- qr(covmat)$rank
p_eff

# Calculate the coefficient of determination (r-squared) for the model:
sstot = sum((PKexample$Concentration-mean(PKexample$Concentration))^2)
sserr = sum((PKexample$Concentration - PKfity)^2)
rsq = 1-sserr/sstot
rsq

#Calculate the Adjusted r^2 for the model fit:
n <- length(PKexample$Concentration)
rsq_adj <- 1 - (1 - rsq) * (n - 1) / (n - p_eff - 1)
rsq_adj

#Raw values of correlation matrix:
cormat

#Plot the correlation matrix (nice colour plot)
#Look for off-diagonal numbers > 0.8 or <-0.8 here. That's bad!
library(corrplot)
corrplot(cormat, method = "color", addCoef.col = "black")
