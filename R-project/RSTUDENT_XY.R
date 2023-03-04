# R-Studentized Residuals Test on Linear Regression

#1) Copy the data set from Excel onto the clipboard (Column: X, then Y, include X and Y as headers)
#2) Make sure the outliers library is installed (Tools -> Install Packages -> "outliers")
#2) Now type in this line and hit Enter (without the # sign)
#For PC:
#DataXY <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#DataXY <- read.table(pipe("pbpaste"), header=TRUE)

attach(DataXY)

#Fit the model using glm (general least-squares mean)
fit1 <- glm(Y~X, data=DataXY, na.action=na.exclude)

#to show residuals and r-studentized residuals:
rs.fit1 <- data.frame(X=X,Y=Y,Residual=residuals(fit1),rStudent=rstudent(fit1))
rs.fit1
