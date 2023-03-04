# Grubb's Test on a Repeated Measurement

#1) Copy the data set from Excel onto the clipboard (Only one column, label X on the top)
#2) Make sure the outliers library is installed (Tools -> Install Packages -> "outliers")
#2) Now type in this line and hit Enter (without the # sign)
#For PC:
#DataX <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#DataX <- read.table(pipe("pbpaste"), header=TRUE)

attach(DataX)

shapiro.test(X) # are the data normally distributed?

shapiro.test(log(X)) # are the data log-normally distributed?
