#This script fits PK data following oral administration with 1-compartment elimination
#to the appropriate model. Written by: David Dubins

#Read in the data (Time, Concentration) for a Subject
#For PC:
#PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

attach(PKexample)

#Plot the log-Concentration vs. Time plot:
plot(Time,log(Concentration),main="Concentration vs. Time: Semi-log Plot", xlab="Time (hours)", ylab = "log(Concentration)")

#An easy way to start is to define the model C(t,1) as a function called Ct1:

Ct1 <- function(A,ka,ke) {
A*(ka/(ke-ka))*(exp(-ka*Time)-exp(-ke*Time))
}

#Now, we never have to write the equation again in the program.

#Define the chisq function as (Yobs-Ymodel)^2 (we are going to minimize this)

chisq <- function(p) sum((Concentration - Ct1(p[1],p[2],p[3]))^2)

#We are using "p" as an array of parameters:
#p[1] = A
#p[2] = ka
#p[3] = ke

#Now guess some values on a new plot.
#Set up variable names g1 to g5 --> easier to change 1st guesses.
#(replace 1, 2, and 3 with your starting guesses)

plot(Time,Concentration,main="Concentration vs. Time", xlab="Time (hours)", ylab = "Concentration (ng/mL)")
g1 <- 25    #first guess for A
g2 <- 1.5    #first guess for ka
g3 <- 0.0488984    #first guess for ke
PKfitx <- PKexample$Time
PKfity <- Ct1(g1,g2,g3)
lines(spline(PKfitx, PKfity))

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2, and 3 with your best guesses)

fit1 <- nlm(chisq, p = c(g1, g2, g3), hessian = TRUE)

#To see the results of the fitting routine:
fit1

#Now superimpose our minimized fit on a brand new plot:
plot(Time,Concentration,main="Concentration vs. Time", xlab="Time (hours)", ylab = "Concentration (ng/mL)")
PKfitx <- PKexample$Time
PKfity <- Ct1(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3])
lines(spline(PKfitx, PKfity))

#Now superimpose our minimized fit on a brand new log plot:
plot(Time,log(Concentration),main="ln(Concentration) vs. Time", xlab="Time (hours)", ylab = "ln(Concentration) (ng/mL)")
PKfitx <- PKexample$Time
PKfitlny <- log(Ct1(fit1$estimate[1],fit1$estimate[2],fit1$estimate[3]))
PKfitlny[1]=NA  #This is to avoid crashing when plotting log(0)
lines(spline(PKfitx, PKfitlny))

#To obtain the approximate standard errors of the parameter estimates:
stderr <- sqrt(diag(2*fit1$minimum/(length(Concentration) - 2) * solve(fit1$hessian)))
stderr

#Calculate the rsq:
sstot = sum((Concentration-mean(Concentration))^2)
sserr = sum((Concentration - PKfity)^2)
rsq = 1-sserr/sstot
rsq

#absorption half-life: (hours)
log(2)/fit1$estimate[2]

#elimination half-life: (hours)
log(2)/fit1$estimate[3]