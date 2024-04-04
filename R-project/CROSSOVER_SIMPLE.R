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
lm.lnaucinf <- lme(-log(AUCinf)~Treatment+Period+Sequence, data=ABEexample1, random=~1|Subject,na.action=na.exclude)
lm.lncmax <- lme(-log(Cmax)~Treatment+Period+Sequence, data=ABEexample1, random=~1|Subject,na.action=na.exclude)

#Analysis for AUCt
summary(lm.lnauct)
intervals(lm.lnauct,0.90,"fixed")
VarCorr(lm.lnauct)

#Analysis for AUCinf
summary(lm.lnaucinf)
intervals(lm.lnaucinf,0.90,"fixed")
VarCorr(lm.lnaucinf)

#Analysis for Cmax
summary(lm.lncmax)
intervals(lm.lncmax,0.90,"fixed")
VarCorr(lm.lncmax)

#roll up the answers
Parameter <- c("lnAUCt","lnAUCinf","lnCmax")
AUCt.stderr <- attr(lm.lnauct$fixDF,"varFixFact")[2,2]
AUCt.intervals <- intervals(lm.lnauct,0.90,"fixed")
AUCt.var <- VarCorr(lm.lnauct)
AUCinf.stderr <- attr(lm.lnaucinf$fixDF,"varFixFact")[2,2]
AUCinf.intervals <- intervals(lm.lnaucinf,0.90,"fixed")
AUCinf.var <- VarCorr(lm.lnaucinf)
Cmax.stderr <- attr(lm.lncmax$fixDF,"varFixFact")[2,2]
Cmax.intervals <- intervals(lm.lncmax,0.90,"fixed")
Cmax.var <- VarCorr(lm.lncmax)

Stderr <- c(AUCt.stderr, AUCinf.stderr, Cmax.stderr)
LowerCI <- c(AUCt.intervals$fixed[2,1], AUCinf.intervals$fixed[2,1],Cmax.intervals$fixed[2,1])
Ratio <- c(AUCt.intervals$fixed[2,2], AUCinf.intervals$fixed[2,2],Cmax.intervals$fixed[2,2])
UpperCI <- c(AUCt.intervals$fixed[2,3], AUCinf.intervals$fixed[2,3],Cmax.intervals$fixed[2,3])

lnResults <- data.frame(Parameter, Stderr, Ratio, LowerCI, UpperCI)

Parameter <- c("AUCt","AUCinf","Cmax")
IntraCV <- c(100*((exp(as.numeric(AUCt.var[2]))-1)^0.5),100*((exp(as.numeric(AUCinf.var[2]))-1)^0.5),100*((exp(as.numeric(Cmax.var[2]))-1)^0.5))
InterCV <- c(100*((exp(as.numeric(AUCt.var[1]))-1)^0.5),100*((exp(as.numeric(AUCinf.var[1]))-1)^0.5),100*((exp(as.numeric(Cmax.var[1]))-1)^0.5))
LowerCI <- 100*exp(LowerCI)
Ratio <- 100*exp(Ratio)
UpperCI <- 100*exp(UpperCI)

Results <- data.frame(Parameter, IntraCV, InterCV, Ratio, LowerCI, UpperCI)

#Now report the results:
lnResults
Results

