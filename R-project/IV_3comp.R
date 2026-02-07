#iv_3comp.R
#This script fits PK data following IV administration with a 2-compartment model.
#Written by: David Dubins
#Date: 04-Feb-26

#Read in the data (Time, Concentration) for a Subject
#For PC:
PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

Dose <- 400 #dose in microlitres

# Define the model you are using to fit the data:
iv_3comp <- function(Vd,L1,L2,L3,k21,k31) {
  (Dose/Vd)*(((k21-L1)*(k31-L1)/((L2-L1)*(L3-L1)))*exp(-L1*PKexample$Time)+((k21-L2)*(k31-L2)/((L1-L2)*(L3-L2)))*exp(-L2*PKexample$Time)+((k21-L3)*(k31-L3)/((L1-L3)*(L2-L3)))*exp(-L3*PKexample$Time))
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((PKexample$Concentration-iv_3comp(p[1],p[2],p[3],p[4],p[5],p[6]))^2)
#We are using "p" as an array of parameters
#p[1] = Vd (in L)
#p[2] = lambda1 (1/h)
#p[3] = lambda2 (1/h)
#p[4] = lambda3 (1/h)
#p[5] = k21 (1/h)
#p[6] = k31 (1/h)

#Plot the original data on a graph:
plot(PKexample$Time,PKexample$Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (ng/mL)")

# Try your hand at fitting your parameters to get in the right ballpark. 
# Myfity has your first guesses at g1, g2, g3, g4, g5, g6, g7 (replace 1, 2, 3, 4, 5, 6, 7)
# with your first guesses)
g1 <- 44.59      #first guess for Vd (L)
g2 <- 6          #first guess for L1
g3 <- 0.18       #first guess for L2
g4 <- 0.01       #first guess for L3
g5 <- 4.54       #first guess for k21
g6 <- 0.11       #first guess for k31
  
# Plot your guess on the plot:
Myfitx <- PKexample$Time
Myfity <- iv_3comp(g1,g2,g3,g4,g5,g6)
lines(spline(Myfitx, Myfity))
#if your first guesses are off the plot, find them:
#plot(Myfitx, Myfity)

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2, g3, g4, g5, g6), hessian = TRUE)
# See the results of your fit:
Myfitout

#Now superimpose our minimized fit on a brand new plot:
plot(PKexample$Time,PKexample$Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (ng/mL)")
PKfitx <- PKexample$Time
PKfity <- iv_3comp(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6])
lines(spline(PKfitx, PKfity))

#Now superimpose our minimized fit on a semi-log plot:
plot(PKexample$Time,log(PKexample$Concentration),main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (ng/mL)")
PKfitx <- PKexample$Time
PKfity_ln <- log(iv_3comp(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6]))
lines(spline(PKfitx, PKfity_ln))

#Calculate the Covariance matrix and the correlation matrix
covmat <- 2 * Myfitout$minimum / (length(PKexample$Concentration) - 2) * solve(Myfitout$hessian)
#Convert Covariance Matrix to the Correlation Matrix.
cormat <- cov2cor(covmat)
param.names <- c("Vd", "L1", "L2", "L3", "k21", "k31")
dimnames(cormat) <- list(param.names, param.names)

# Calculate the standard errors of parameter estimates
semat <- sqrt(diag(covmat))
names(semat) <- c("Vd (L)", "L1 (1/h)", "L2 (1/h)", "L3 (1/h)", "k21 (1/h)", "k31 (1/h)")

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
