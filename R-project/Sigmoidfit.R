#This script fits a sigmoid model to data (X values, Y values)
#Written by: David Dubins

#For PC:
#Mydata <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#Mydata <- read.table(pipe("pbpaste"), header=TRUE)
attach(Mydata)
# Define the model you are using to fit the data (in this case, sigmoid):
Sigfn <- function(A,B,C,D) {
(D+(A/(1+exp(-B*X+C))))
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((Y-Sigfn(p[1],p[2],p[3],p[4]))^2)
# Now plot your data:
plot(X,Y)

# Try your hand at fitting your parameters to get in the right ballpark. 
# Myfity has your first guesses at g1, g2, g3, and g4 (replace 1, 2, 3, 4
# with your first guesses)

g1 <- 9.5 # first guess for parameter A
g2 <- 0.5 # first guess for parameter B
g3 <- 10 # first guess for parameter C
g4 <- 3 # first guess for parameter D

Myfitx <- Mydata$X
Myfity <- Sigfn(g1,g2,g3,g4)

# Plot your guess on the plot:
lines(smooth.spline(Myfitx, Myfity, spar = 0))   #spar is smoothing param (0.0: no smoothing, 1.0:too far)

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2, g3, g4), hessian = TRUE)
# See the results of your fit:
Myfitout

# To obtain the approximate standard errors of the parameter estimates:
se <- sqrt(diag(2*Myfitout$minimum/(length(Y) - 2) * solve(Myfitout$hessian)))
se

# Now superimpose our minimized fit on a brand new plot:
plot(X,Y)
Myfitx <- Mydata$X
Myfity <- Sigfn(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4])
lines(smooth.spline(Myfitx, Myfity, spar=0))

#Now plot the residuals to see if they look reasonable:
Resid <- Mydata$Y - Sigfn(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4])
plot(X,Resid,main="Residuals vs. X", xlab="Time (hours)", ylab = "Yobserved - Yfit")
abline(h = 0) # plot a line at y=0

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Y-mean(Y))^2)
sserr = sum((Y - Myfity)^2)
rsq = 1-sserr/sstot
rsq

# Now calculate the midpoint of the transition, (MidX, MidY) (=C/B)
MidX <- Myfitout$estimate[3]/Myfitout$estimate[2]
MidY <-  (Myfitout$estimate[1]/(1+exp(-Myfitout$estimate[2]*MidX+Myfitout$estimate[3])))+Myfitout$estimate[4]


#Report the midpoint of the transition on the X-axis and the value of Y at this point:
MidX
MidY

#Report the approximate standard error of X (using partial derivative weighted errors)
se_MidX <- sqrt((se[3]/Myfitout$estimate[2])^2+(se[2]*Myfitout$estimate[3]/(Myfitout$estimate[2]^2))^2)
se_MidX
