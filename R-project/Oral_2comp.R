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
A*ka*( ((k21-ka)*exp(-ka*Time)/((L1-ka)*(L2-ka)))
      +((k21-L1)*exp(-L1*Time)/((ka-L1)*(L2-L1)))
      +((k21-L2)*exp(-L2*Time)/((ka-L2)*(L1-L2))) )
}

#Now, we never have to write the equation again in the program.

#Define the function as (Yobs-Ymodel)^2 (we are going to minimize this)

chisq2 <- function(p) sum((Concentration - Ct1(p[1],p[2],p[3],p[4],p[5]))^2)

#We are using "p" as an array of parameters:
#p[1] = A
#p[2] = ka
#p[3] = k21
#p[4] = L1
#p[5] = L2

#Unfortunately, most fitting programs need a starting point. You can't
#just say "fit this equation" and expect it to do it, you have to start
#with values in the right ball park.

#Guess some values on a new plot 
#Set up variable names g1 to g5 --> easier to change 1st guesses
#(replace 1, 2, 3, 4, 5 with your starting guesses)
plot(Time,Concentration,main="Concentration vs. Time", xlab="Time (hours)", ylab = "Concentration (ng/mL)")
g1 <- 1  #first guess for A
g2 <- 2  #first guess for ka
g3 <- 3  #first guess for k21
g4 <- 4  #first guess for L1
g5 <- 5  #first guess for L2
PKfitx <- PKexample$Time
PKfity <- Ct1(g1,g2,g3,g4,g5)
lines(spline(PKfitx, PKfity))

#Do do the actual fit with some trivial starting values. This is done
#by minimizing the function we defined using the starting values we just tried:
#(replace 1, 2, 3, 4, 5 with your best guesses)

fit2 <- nlm(chisq2, p = c(1, 2, 3, 4, 5), hessian = TRUE)

#To see the results of the fitting routine:
fit2

#To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*fit2$minimum/(length(Concentration) - 2) * solve(fit2$hessian)))

#Now superimpose our minimized fit on a brand new plot:
plot(Time,Concentration,main="Concentration vs. Time", xlab="Time (hours)", ylab = "Concentration (ng/mL)")
PKfitx <- PKexample$Time
PKfity <- Ct1(fit2$estimate[1],fit2$estimate[2],fit2$estimate[3],fit2$estimate[4],fit2$estimate[5])
lines(spline(PKfitx, PKfity))

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Fwhatever-mean(Fwhatever))^2)
sserr = sum((Fwhatever - Myfity)^2)
rsq = 1-sserr/sstot
rsq