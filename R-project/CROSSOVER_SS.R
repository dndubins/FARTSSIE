#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental_PK_SS" (including headers), Control-C
#3) Now type in this line and hit Enter:

#For PC:
#SSexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#SSexample1 <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R:
library(nlme)
factors       <- c("Subject", "Sequence", "Period", "Treatment")
ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for lme()
str(ABEexample1) #check the data

lm.lnAUCtau <- lme(-log(AUCtau)~Treatment+Period+Sequence, data=SSexample1, random=~1|Subject,na.action=na.exclude)
lm.lnCssmax <- lme(-log(Cssmax)~Treatment+Period+Sequence, data=SSexample1, random=~1|Subject,na.action=na.exclude)
lm.lnCssmin <- lme(-log(Cssmin)~Treatment+Period+Sequence, data=SSexample1, random=~1|Subject,na.action=na.exclude)
lm.Fluctuation <- lme(Fluctuation~Treatment+Period+Sequence, data=SSexample1, random=~1|Subject,na.action=na.exclude)

#NOTE: If you don't get intervals, use fixed intervals(lm.varname,0.90,"fixed")
#Analysis for AUCtau
summary(lm.lnAUCtau)
intervals(lm.lnAUCtau,0.90)
VarCorr(lm.lnAUCtau)

#Analysis for Cssmax
summary(lm.lnCssmax)
intervals(lm.lnCssmax,0.90)
VarCorr(lm.lnCssmax)

#Analysis for Cssmin
summary(lm.lnCssmin)
intervals(lm.lnCssmin,0.90)
VarCorr(lm.lnCssmin)

#Analysis for Fluctuation
summary(lm.Fluctuation)
intervals(lm.Fluctuation,0.90)
VarCorr(lm.Fluctuation)
