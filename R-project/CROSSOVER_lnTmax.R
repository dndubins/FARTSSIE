#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers), Control-C
#3) Now type in this line and hit Enter (without the # sign)
#ABEexample1 <- read.table(file="clipboard", header=TRUE)

#Now copy and paste the following into R: (or select "Edit", then "Run all")

library(nlme)
factors       <- c("Subject", "Sequence", "Period", "Treatment")
ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for lme()
str(ABEexample1) #check the data

lm.lnTmax <- lme(-log(Tmax)~Treatment+Period+Sequence, data=ABEexample1, random=~1|Subject,na.action=na.exclude)

#Analysis for Tmax
summary(lm.lnTmax)
intervals(lm.lnTmax,0.90,"fixed")
VarCorr(lm.lnTmax)

#roll up the answers
Parameter <- c("lnTmax")
Tmax.stderr <- attr(lm.lnTmax$fixDF,"varFixFact")[2,2]
Tmax.intervals <- intervals(lm.lnTmax,0.90,"fixed")
Tmax.var <- VarCorr(lm.lnTmax)

Stderr <- c(Tmax.stderr)
LowerCI <- c(Tmax.intervals$fixed[2,1])
Ratio <- c(Tmax.intervals$fixed[2,2])
UpperCI <- c(Tmax.intervals$fixed[2,3])

lnResults <- data.frame(Parameter, Stderr, Ratio, LowerCI, UpperCI)

Parameter <- c("Tmax")
IntraCV <- c(100*((exp(as.numeric(Tmax.var[2]))-1)^0.5))
InterCV <- c(100*((exp(as.numeric(Tmax.var[1]))-1)^0.5))
LowerCI <- 100*exp(LowerCI)
Ratio <- 100*exp(Ratio)
UpperCI <- 100*exp(UpperCI)

Results <- data.frame(Parameter, IntraCV, InterCV, Ratio, LowerCI, UpperCI)

#Now report the results:
lnResults
Results
