#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers), Control-C
#3) Now type in this line and hit Enter (without the # sign)
#For PC:
#ABEexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R: (or select "Edit", then "Run all")

library(nlme)
factors       <- c("Subject", "Sequence", "Period", "Treatment")
ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for lme()
str(ABEexample1) #check the data

lm.lnauct <- lme(-log(AUCt)~Treatment+Period+Sequence, data=ABEexample1, random=~1|Subject,na.action=na.exclude)

#Analysis for AUCt
summary(lm.lnauct)
intervals(lm.lnauct,0.9,"fixed")
VarCorr(lm.lnauct)

#roll up the answers
Parameter <- c("lnAUCt")
AUCt.stderr <- attr(lm.lnauct$fixDF,"varFixFact")[2,2]
AUCt.intervals <- intervals(lm.lnauct,0.9,"fixed")
AUCt.var <- VarCorr(lm.lnauct)

Stderr <- c(AUCt.stderr)
LowerCI <- c(AUCt.intervals$fixed[2,1])
Ratio <- c(AUCt.intervals$fixed[2,2])
UpperCI <- c(AUCt.intervals$fixed[2,3])

lnResults <- data.frame(Parameter, Stderr, Ratio, LowerCI, UpperCI)

Parameter <- c("AUCt")
IntraCV <- c(100*((exp(as.numeric(AUCt.var[2]))-1)^0.5))
InterCV <- c(100*((exp(as.numeric(AUCt.var[1]))-1)^0.5))
LowerCI <- 100*exp(LowerCI)
Ratio <- 100*exp(Ratio)
UpperCI <- 100*exp(UpperCI)
n <- lm.lnauct$dims$N/2
tau1 <- sqrt(n)*log(95/80)/(sqrt(2*as.numeric(AUCt.var[2])))
tau2 <- sqrt(n)*log(95/125)/(sqrt(2*as.numeric(AUCt.var[2])))
Power <- 100*(pt(tau1-qt(0.95,n-2),df=n-2) - pt(tau2+qt(0.95,n-2),df=n-2))
#Note that the canadian guidance doesn't adjust for the non-centrality parameter for the calculation of power.
#Also note that Method C (Povtin 2008) calculates the power based on T/R = 0.95.
Results <- data.frame(Parameter, IntraCV, InterCV, Ratio, LowerCI, UpperCI, Power)

#Now report the results:
lnResults
Results

