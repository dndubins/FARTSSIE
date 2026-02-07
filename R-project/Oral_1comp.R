#Oral_1comp.R
#This script fits PK data following oral administration with 1-compartment elimination
#Written by: David Dubins
#Date: 05-Feb-26

#Read in the data (Time, Concentration) for a Subject
#For PC:
#PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

#Save a bit of referencing work:
#(otherwise we would need to write PKexample$Time and PKexample$Concentration everywhere)
Time <- PKexample$Time 
Concentration <- PKexample$Concentration

#Have a look at the data to make sure it was imported correctly:
plot(Time,Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (µg/mL)")

#Define the model C(t,1) as a function called Ct1:
Ct1 <- function(A,ka,ke) {
A*(ka/(ke-ka))*(exp(-ka*Time)-exp(-ke*Time))
}

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)
chisq <- function(p) sum((Concentration - Ct1(p[1],p[2],p[3]))^2)
#We are using "p" as an array of parameters:
#p[1] = A
#p[2] = ka
#p[3] = ke

#Now guess some values on a new plot.
#Set up variable names g1 to g3 --> easier to change 1st guesses.
#(replace 1, 2, and 3 with your starting guesses)

plot(Time,Concentration,main = "Concentration vs. Time", xlab = "Time (h)", ylab = "Concentration (µg/mL)")
g1 <- 1    #first guess for A
g2 <- 1    #first guess for ka
g3 <- 1    #first guess for ke
PKfitx <- Time
PKfity <- Ct1(g1,g2,g3)
lines(smooth.spline(PKfitx, PKfity, spar = 0)) #spar is smoothing param (0.0: no smoothing, 1.0:too far)

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2, and 3 with your best guesses)
fit1 <- nlm(chisq, p = c(g1, g2, g3), hessian = TRUE)

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand-new plot:
plot(Time,Concentration,main="Concentration vs. Time", xlab="Time (h)", ylab = "Concentration (ng/mL)")
PKfitx <- Time
PKfity <- Ct1(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3])
lines(smooth.spline(PKfitx, PKfity, spar=0))

#Now superimpose our minimized fit on a brand-new log plot:
plot(Time,log(Concentration),main="ln(Concentration) vs. Time", xlab="Time (h)", ylab = "ln(Concentration) (ln(ng/mL))")
PKfitx <- Time
PKfitlny <- log(Ct1(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3]))
PKfitlny[1]=0  #This is to avoid crashing when plotting log(0)
lines(smooth.spline(PKfitx, PKfitlny, spar = 0))

#Now plot the residuals to see if they look reasonable:
Resid <- PKexample$Concentration-Ct1(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3])
plot(Time,Resid,main="Residuals vs. Time", xlab="Time (h)", ylab = "Cobserved - Cfit")
abline(h = 0)  #plot a line at y=0

#Solve for the correlation matrix:
cov.mat <- 2 * fit1$minimum / (length(Concentration) - 2) * solve(fit1$hessian)
rownames(cov.mat) <- colnames(cov.mat) <- c("A","ka","ke")
cor.mat <- cov2cor(cov.mat)
cor.mat

# Calculate the approximate standard errors of the parameter estimates:
stderr <- sqrt(diag(cov.mat))
stderr

#Calculate the coefficient of determination (r-squared) for the model:
sstot = sum((Concentration-mean(Concentration))^2)
sserr = sum((Concentration - PKfity)^2)
rsq = 1-sserr/sstot
rsq

#absorption half-life: (hours)
log(2)/fit1$estimate[2]

#elimination half-life: (hours)
log(2)/fit1$estimate[3]

#Plot the correlation matrix (nice colour plot)
#Look for off-diagonal numbers > 0.8 or <-0.8 here. That's bad!
library(corrplot)
corrplot(cor.mat, method = "color", addCoef.col = "black", number.digits = 3)
