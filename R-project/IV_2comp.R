#This script fits PK data following IV administration with 2-compartment elimination
#to the appropriate model. Written by: David Dubins

#Read in the data (Time, Concentration) for a Subject
#For PC:
#PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

attach(PKexample)

#Define the function as (Yobs-Ymodel)^2 (we are going to minimize this)

fn <- function(p) sum((Concentration - (p[1]*(p[2]-p[3])*(exp(-p[3]*Time))-p[1]*(p[2]-p[4])*exp(-p[4]*Time)))^2)

#The model itself is Concentration=A(exp(-ka*Time)-exp(-ke*Time)
#but we are re-arranging the expression to:
#(Concentration - A(exp(-ka*time)-exp(-ke*Time)))^2  <--- we're minimizing this
#and we are using "p" as an array of parameters instead of A, ka, and ke:
#p[1] = A
#p[2] = k21
#p[3] = lambda1
#p[4] = lambda2

#Plot the original data on a graph:
plot(Time,Concentration)

#Unfortunately, most fitting programs need a starting point. You can't
#just say "fit this equation" and expect it to do it, you have to start
#with values in the right ball park.

#Now try some initial values and see what they look like.

PKfitx <- PKexample$Time
PKfity <- -20*(0.05-1)*(exp(-1*Time))+20*(0.05-0.03)*exp(-0.03*Time)
lines(spline(PKfitx, PKfity))

#Do do the actual fit with the starting values we found. This is done
#by minimizing the function we defined using the starting values we just tried:
PKfitout <- nlm(fn, p = c(-20, 0.05, 1, 0.03), hessian = TRUE)

#To see the results of the fitting routine:
PKfitout

#To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*PKfitout$minimum/(length(Concentration) - 2) * solve(PKfitout$hessian)))

#Now superimpose our minimized fit on a brand new plot:
plot(Time,Concentration)
PKfitx <- PKexample$Time
PKfity <- PKfitout$estimate[1]*(PKfitout$estimate[2]-PKfitout$estimate[3])*(exp(-PKfitout$estimate[3]*Time))-PKfitout$estimate[1]*(PKfitout$estimate[2]-PKfitout$estimate[4])*exp(-PKfitout$estimate[4]*Time)

lines(spline(PKfitx, PKfity))

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Fwhatever-mean(Fwhatever))^2)
sserr = sum((Fwhatever - Myfity)^2)
rsq = 1-sserr/sstot
rsq

