#This script fits PK data following IV administration with 3-compartment elimination
#to the appropriate model. Written by: David Dubins

#Read in the data (Time, Concentration) for a Subject
#For PC:
PKexample <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#PKexample <- read.table(pipe("pbpaste"), header=TRUE)

attach(PKexample)
Dose <- 400 #dose in microlitres

#Define the function as (Yobs-Ymodel)^2 (we are going to minimize this)

# Define the model you are using to fit the data (in this case, sigmoid):
fn <- function(L1,L2,L3,k21,k31,Vd) {
  (Dose/Vd)*(((k21-L1)*(k31-L1)/((L2-L1)*(L3-L1)))*exp(-L1*Time)+((k21-L2)*(k31-L2)/((L1-L2)*(L3-L2)))*exp(-L2*Time)+((k21-L3)*(k31-L3)/((L1-L3)*(L2-L3)))*exp(-L3*Time))
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((Concentration-fn(p[1],p[2],p[3],p[4],p[5],p[6]))^2)

#The model itself is Concentration=something very long
#We are using "p" as an array of parameters instead of A, ka, and ke:
#p[1] = L1
#p[2] = L2
#p[3] = L3
#p[4] = k21
#p[5] = k31
#p[6] = Vd

#Plot the original data on a graph:
plot(Time,Concentration)

# Try your hand at fitting your parameters to get in the right ballpark. 
# Myfity has your first guesses at g1, g2, g3, g4, g5, g6, g7 (replace 1, 2, 3, 4, 5, 6, 7)
# with your first guesses)

g1 <- 6     #first guess for L1
g2 <- 0.18  #first guess for L2
g3 <- 0.01  #first guess for L3
g4 <- 4.54   #first guess for k21
g5 <- 0.11   #first guess for k31
g6 <- 44.59    #first guess for Vd

Myfitx <- PKexample$Time
Myfity <- fn(g1,g2,g3,g4,g5,g6)

# Plot your guess on the plot:
lines(spline(Myfitx, Myfity))

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2, g3, g4, g5, g6), hessian = TRUE)
# See the results of your fit:
Myfitout


#To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*Myfitout$minimum/(length(Concentration) - 2) * solve(Myfitout$hessian)))

#Now superimpose our minimized fit on a brand new plot:
plot(Time,Concentration)
PKfitx <- PKexample$Time
PKfity <- fn(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6])
lines(spline(PKfitx, PKfity))

#Now superimpose our minimized fit on a semi-log plot:
plot(Time,log(Concentration))
PKfitx <- PKexample$Time
PKfity_ln <- log(fn(Myfitout$estimate[1],Myfitout$estimate[2],Myfitout$estimate[3],Myfitout$estimate[4],Myfitout$estimate[5],Myfitout$estimate[6]))
lines(spline(PKfitx, PKfity_ln))

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Concentration-mean(Concentration))^2)
sserr = sum((Concentration - PKfity)^2)
rsq = 1-sserr/sstot
rsq

