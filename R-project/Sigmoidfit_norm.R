#Sigmoidfit_norm.R
#This script fits a sigmoid model to normalized data (X values, Y values (0-1))
#Written by: David Dubins
#Date: 06-Feb-26

#For PC:
#Mydata <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#Mydata <- read.table(pipe("pbpaste"), header=TRUE)

#Save a bit of referencing work:
#(otherwise we would need to write Mydata$X and Mydata$Ynorm everywhere)
X <- Mydata$X
Ynorm <- Mydata$Ynorm

# Define the model you are using to fit the data (in this case, sigmoid):
Sigfn <- function(B,C) {
  (1/(1+exp(-B*(X-C))))
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((Ynorm-Sigfn(p[1],p[2]))^2)
# Now plot your data:
plot(X,Ynorm)

# Try your hand at fitting your parameters to get in the right ballpark. 
# Myfity has your first guesses at g1, g2

g1 <- 0.4 # first guess for parameter B
g2 <- 21 # first guess for parameter C

Myfitx <- X
Myfity <- Sigfn(g1,g2)

# Plot your guess on the plot:
lines(smooth.spline(Myfitx, Myfity, spar = 0))   #spar is smoothing param (0.0: no smoothing, 1.0:too far)

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2), hessian = TRUE)
# See the results of your fit:
Myfitout

# Now superimpose our minimized fit on a brand new plot:
plot(X,Ynorm)
Myfitx <- Mydata$X
Myfity <- Sigfn(Myfitout$estimate[1],Myfitout$estimate[2])
lines(smooth.spline(Myfitx, Myfity, spar=0))

#Now plot the residuals to see if they look reasonable:
Resid <-Ynorm-Sigfn(Myfitout$estimate[1],Myfitout$estimate[2])
plot(X,Resid,main="Residuals vs. Time", xlab="X", ylab = "Yobserved - Yfit")
abline(h = 0)  #plot a line at y=0

#Solve for the correlation matrix:
cov.mat <- 2 * Myfitout$minimum / (length(Ynorm) - 2) * solve(Myfitout$hessian)
rownames(cov.mat) <- colnames(cov.mat) <- c("B","C")
cor.mat <- cov2cor(cov.mat)
cor.mat

# Calculate the approximate standard errors of the parameter estimates:
stderr <- sqrt(diag(cov.mat))
stderr

#Calculate the coefficient of determination (r-squared) for the model:
sstot = sum((Ynorm-mean(Ynorm))^2)
sserr = sum((Ynorm - Myfity)^2)
rsq = 1-sserr/sstot
rsq

# Now calculate the midpoint of the transition, (MidX, MidY)
MidX <- Myfitout$estimate[2]
MidY <- (1/(1+exp(-Myfitout$estimate[1]*(MidX-Myfitout$estimate[2]))))

#Report the midpoint of the transition on the X-axis and the value of Y at this point:
MidX
MidY

#Plot the correlation matrix (nice colour plot)
#Look for off-diagonal numbers > 0.8 or <-0.8 here. That's bad!
library(corrplot)
corrplot(cor.mat, method = "color", addCoef.col = "black", number.digits = 3)
