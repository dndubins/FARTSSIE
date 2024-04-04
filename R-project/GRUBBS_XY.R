# Grubb's Test on Linear Regression (residuals)
# Source: https://sites.chem.utoronto.ca/chemistry/coursenotes/analsci/stats/Outliers.html

#1) Copy the data set from Excel onto the clipboard (Column: X, then Y, include X and Y as headers)
#2) Make sure the outliers library is installed (Tools -> Install Packages -> "outliers")
#2) Now type in this line and hit Enter (without the # sign)
#For PC:
#DataXY <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#DataXY <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R: (change variable as necessary)

library(outliers)

attach(DataXY)  #attach the data so we can use header names as variables

# Define whatever model you are using to fit the data (in this case, linear):
fn <- function(A,B) {
  (A+B*X)
}

# Define the chisq function as sum(Yobs - Ymodel)^2:
chisq <- function(p) sum((Y-fn(p[1],p[2]))^2)

# Now plot your data:
plot(X,Y)

# Try some starting guesses
g1 <- 2 #here, g1 is intercept
g2 <- 5 #here, g2 is slope
# Calculate and plot fit based on starting guesses and plot them (repeat as necessary
Myfitx <- X
Myfity <- fn(g1,g2)
lines(spline(Myfitx, Myfity))

# Ready for fitting the model with your guesses:
Myfitout <- nlm(chisq, p = c(g1, g2), hessian = TRUE)
# See the results of your fit:
Myfitout

# To obtain the approximate standard errors of the parameter estimates:
sqrt(diag(2*Myfitout$minimum/(length(Y) - 2) * solve(Myfitout$hessian)))

# Now superimpose our minimized fit on a brand new plot:
plot(X,Y)
Myfitx <- X
Myfity <- fn(Myfitout$estimate[1],Myfitout$estimate[2])
lines(spline(Myfitx, Myfity))

# You can also calculate the coefficient of determination (r-squared) for your model:
sstot = sum((Y-mean(Y))^2)
sserr = sum((Y - Myfity)^2)
rsq = 1-sserr/sstot
rsq

# Now let's plot the residuals (Yobs-Ymodel)
resid = Y-fn(Myfitout$estimate[1],Myfitout$estimate[2])
plot(resid)
# Run the grubb's test on the residuals
grubbs.test(resid, type = 10, opposite = FALSE, two.sided = TRUE)
# "type" is an integer value indicating test variant. 
# 10 is a test for one outlier (side is detected automatically and 
# can be reversed by opposite parameter).
# 11 is a test for two outliers on opposite tails,
# 20 is test for two outliers in one tail. 
