# Grubb's Test on a Repeated Measurement

#1) Copy the data set from Excel onto the clipboard (Only one column, label X on the top)
#2) Make sure the outliers library is installed (Tools -> Install Packages -> "outliers")
#2) Now type in this line and hit Enter (without the # sign)
#For PC:
#DataX <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#DataX <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R: (change variable as necessary)

library(outliers)

attach(DataX)  #attach the data so we can use header names as variables

# Plot the data (have a look!)
plot(X)

# Run the grubb's test on the residuals
grubbs.test(X, type = 10, opposite = FALSE, two.sided = TRUE)
# "type" is an integer value indicating test variant. 
# 10 is a test for one outlier (side is detected automatically and 
# can be reversed by opposite parameter).
# 11 is a test for two outliers on opposite tails,
# 20 is test for two outliers in one tail. 
