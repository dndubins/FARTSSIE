#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers, from "Subject" to "AUCRatio"), Control-C
#3) Now type in this line and hit Enter (without the # sign)
#For PC:
#ABEexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R: (or select "Edit", then "Run all")
attach(ABEexample1)
library(nlme)
lm.lnauct <- lme(-log(AUCt)~Treatment, data=ABEexample1,random=~1|Subject,na.action=na.exclude)
lm.lnaucinf <- lme(-log(AUCinf)~Treatment, data=ABEexample1,random=~1|Subject,na.action=na.exclude)
lm.lncmax <- lme(-log(Cmax)~Treatment, data=ABEexample1,random=~1|Subject,na.action=na.exclude)

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

ANOVA_lnResults <- data.frame(Parameter, Stderr, Ratio, LowerCI, UpperCI)

Parameter <- c("AUCt","AUCinf","Cmax")
InterCV <- c(100*((exp(as.numeric(AUCt.var[1]))-1)^0.5),100*((exp(as.numeric(AUCinf.var[1]))-1)^0.5),100*((exp(as.numeric(Cmax.var[1]))-1)^0.5))
LowerCI <- 100*exp(LowerCI)
Ratio <- 100*exp(Ratio)
UpperCI <- 100*exp(UpperCI)
pVal_TRT <- c(summary(lm.lnauct)$tTable[2,5],summary(lm.lnaucinf)$tTable[2,5],summary(lm.lncmax)$tTable[2,5])

#Calculate Power
n_auct <- lm.lnauct$dims$N
tau1_auct <- sqrt(n_auct)*log(Ratio[1]/80)/(sqrt(2*as.numeric(AUCt.var[1])))
tau2_auct <- sqrt(n_auct)*log(Ratio[1]/125)/(sqrt(2*as.numeric(AUCt.var[1])))
Power_auct <- 100*(pt(tau1_auct-qt(0.95,2*n_auct-2),df=2*n_auct-2) - pt(tau2_auct+qt(0.95,2*n_auct-2),df=2*n_auct-2))

n_aucinf <- lm.lnaucinf$dims$N
tau1_aucinf <- sqrt(n_aucinf)*log(Ratio[2]/80)/(sqrt(2*as.numeric(AUCinf.var[1])))
tau2_aucinf <- sqrt(n_aucinf)*log(Ratio[2]/125)/(sqrt(2*as.numeric(AUCinf.var[1])))
Power_aucinf <- 100*(pt(tau1_aucinf-qt(0.95,2*n_aucinf-2),df=2*n_aucinf-2) - pt(tau2_aucinf+qt(0.95,2*n_aucinf-2),df=2*n_aucinf-2))

n_cmax <- lm.lncmax$dims$N
tau1_cmax <- sqrt(n_cmax)*log(Ratio[3]/80)/(sqrt(2*as.numeric(Cmax.var[1])))
tau2_cmax <- sqrt(n_cmax)*log(Ratio[3]/125)/(sqrt(2*as.numeric(Cmax.var[1])))
Power_cmax <- 100*(pt(tau1_cmax-qt(0.95,2*n_cmax-2),df=2*n_cmax-2) - pt(tau2_cmax+qt(0.95,2*n_cmax-2),df=2*n_cmax-2))

Power <- c(Power_auct, Power_aucinf, Power_cmax)

ANOVA_Results <- data.frame(Parameter, InterCV, Ratio, LowerCI, UpperCI, pVal_TRT, Power)

#Calculate Descriptive Statistics for Treatment A:
Statistic <- c("n","Mean","Median","Min","Max","STDEV","%CV")
AUCt_A <- c(length(na.omit(AUCt[Treatment == "A"])),mean(AUCt[Treatment == "A"],na.rm=TRUE),median(AUCt[Treatment == "A"],na.rm=TRUE),min(AUCt[Treatment == "A"],na.rm=TRUE),max(AUCt[Treatment == "A"],na.rm=TRUE),sd(AUCt[Treatment == "A"],na.rm=TRUE),100*sd(AUCt[Treatment == "A"],na.rm=TRUE)/mean(AUCt[Treatment == "A"],na.rm=TRUE))
AUCinf_A <- c(length(na.omit(AUCinf[Treatment == "A"])),mean(AUCinf[Treatment == "A"],na.rm=TRUE),median(AUCinf[Treatment == "A"],na.rm=TRUE),min(AUCinf[Treatment == "A"],na.rm=TRUE),max(AUCinf[Treatment == "A"],na.rm=TRUE),sd(AUCinf[Treatment == "A"],na.rm=TRUE),100*sd(AUCinf[Treatment == "A"],na.rm=TRUE)/mean(AUCinf[Treatment == "A"],na.rm=TRUE))
Cmax_A <- c(length(na.omit(Cmax[Treatment == "A"])),mean(Cmax[Treatment == "A"],na.rm=TRUE),median(Cmax[Treatment == "A"],na.rm=TRUE),min(Cmax[Treatment == "A"],na.rm=TRUE),max(Cmax[Treatment == "A"],na.rm=TRUE),sd(Cmax[Treatment == "A"],na.rm=TRUE),100*sd(Cmax[Treatment == "A"],na.rm=TRUE)/mean(Cmax[Treatment == "A"],na.rm=TRUE))
Tmax_A <- c(length(na.omit(Tmax[Treatment == "A"])),mean(Tmax[Treatment == "A"],na.rm=TRUE),median(Tmax[Treatment == "A"],na.rm=TRUE),min(Tmax[Treatment == "A"],na.rm=TRUE),max(Tmax[Treatment == "A"],na.rm=TRUE),sd(Tmax[Treatment == "A"],na.rm=TRUE),100*sd(Tmax[Treatment == "A"],na.rm=TRUE)/mean(Tmax[Treatment == "A"],na.rm=TRUE))
Lambda_A <- c(length(na.omit(Lambda[Treatment == "A"])),mean(Lambda[Treatment == "A"],na.rm=TRUE),median(Lambda[Treatment == "A"],na.rm=TRUE),min(Lambda[Treatment == "A"],na.rm=TRUE),max(Lambda[Treatment == "A"],na.rm=TRUE),sd(Lambda[Treatment == "A"],na.rm=TRUE),100*sd(Lambda[Treatment == "A"],na.rm=TRUE)/mean(Lambda[Treatment == "A"],na.rm=TRUE))
Thalf_A <- c(length(na.omit(T1.2[Treatment == "A"])),mean(T1.2[Treatment == "A"],na.rm=TRUE),median(T1.2[Treatment == "A"],na.rm=TRUE),min(T1.2[Treatment == "A"],na.rm=TRUE),max(T1.2[Treatment == "A"],na.rm=TRUE),sd(T1.2[Treatment == "A"],na.rm=TRUE),100*sd(T1.2[Treatment == "A"],na.rm=TRUE)/mean(T1.2[Treatment == "A"],na.rm=TRUE))
AUCRatio_A <- c(length(na.omit(AUCRatio[Treatment == "A"])),mean(AUCRatio[Treatment == "A"],na.rm=TRUE),median(AUCRatio[Treatment == "A"],na.rm=TRUE),min(AUCRatio[Treatment == "A"],na.rm=TRUE),max(AUCRatio[Treatment == "A"],na.rm=TRUE),sd(AUCRatio[Treatment == "A"],na.rm=TRUE),100*sd(AUCRatio[Treatment == "A"],na.rm=TRUE)/mean(AUCRatio[Treatment == "A"],na.rm=TRUE))

Results_A <- data.frame(Statistic, AUCt_A, AUCinf_A, Cmax_A, Tmax_A, Lambda_A, Thalf_A, AUCRatio_A)

#Calculate Descriptive Statistics for Treatment B:
AUCt_B <- c(length(na.omit(AUCt[Treatment == "B"])),mean(AUCt[Treatment == "B"],na.rm=TRUE),median(AUCt[Treatment == "B"],na.rm=TRUE),min(AUCt[Treatment == "B"],na.rm=TRUE),max(AUCt[Treatment == "B"],na.rm=TRUE),sd(AUCt[Treatment == "B"],na.rm=TRUE),100*sd(AUCt[Treatment == "B"],na.rm=TRUE)/mean(AUCt[Treatment == "B"],na.rm=TRUE))
AUCinf_B <- c(length(na.omit(AUCinf[Treatment == "B"])),mean(AUCinf[Treatment == "B"],na.rm=TRUE),median(AUCinf[Treatment == "B"],na.rm=TRUE),min(AUCinf[Treatment == "B"],na.rm=TRUE),max(AUCinf[Treatment == "B"],na.rm=TRUE),sd(AUCinf[Treatment == "B"],na.rm=TRUE),100*sd(AUCinf[Treatment == "B"],na.rm=TRUE)/mean(AUCinf[Treatment == "B"],na.rm=TRUE))
Cmax_B <- c(length(na.omit(Cmax[Treatment == "B"])),mean(Cmax[Treatment == "B"],na.rm=TRUE),median(Cmax[Treatment == "B"],na.rm=TRUE),min(Cmax[Treatment == "B"],na.rm=TRUE),max(Cmax[Treatment == "B"],na.rm=TRUE),sd(Cmax[Treatment == "B"],na.rm=TRUE),100*sd(Cmax[Treatment == "B"],na.rm=TRUE)/mean(Cmax[Treatment == "B"],na.rm=TRUE))
Tmax_B <- c(length(na.omit(Tmax[Treatment == "B"])),mean(Tmax[Treatment == "B"],na.rm=TRUE),median(Tmax[Treatment == "B"],na.rm=TRUE),min(Tmax[Treatment == "B"],na.rm=TRUE),max(Tmax[Treatment == "B"],na.rm=TRUE),sd(Tmax[Treatment == "B"],na.rm=TRUE),100*sd(Tmax[Treatment == "B"],na.rm=TRUE)/mean(Tmax[Treatment == "B"],na.rm=TRUE))
Lambda_B <- c(length(na.omit(Lambda[Treatment == "B"])),mean(Lambda[Treatment == "B"],na.rm=TRUE),median(Lambda[Treatment == "B"],na.rm=TRUE),min(Lambda[Treatment == "B"],na.rm=TRUE),max(Lambda[Treatment == "B"],na.rm=TRUE),sd(Lambda[Treatment == "B"],na.rm=TRUE),100*sd(Lambda[Treatment == "B"],na.rm=TRUE)/mean(Lambda[Treatment == "B"],na.rm=TRUE))
Thalf_B <- c(length(na.omit(T1.2[Treatment == "B"])),mean(T1.2[Treatment == "B"],na.rm=TRUE),median(T1.2[Treatment == "B"],na.rm=TRUE),min(T1.2[Treatment == "B"],na.rm=TRUE),max(T1.2[Treatment == "B"],na.rm=TRUE),sd(T1.2[Treatment == "B"],na.rm=TRUE),100*sd(T1.2[Treatment == "B"],na.rm=TRUE)/mean(T1.2[Treatment == "B"],na.rm=TRUE))
AUCRatio_B <- c(length(na.omit(AUCRatio[Treatment == "B"])),mean(AUCRatio[Treatment == "B"],na.rm=TRUE),median(AUCRatio[Treatment == "B"],na.rm=TRUE),min(AUCRatio[Treatment == "B"],na.rm=TRUE),max(AUCRatio[Treatment == "B"],na.rm=TRUE),sd(AUCRatio[Treatment == "B"],na.rm=TRUE),100*sd(AUCRatio[Treatment == "B"],na.rm=TRUE)/mean(AUCRatio[Treatment == "B"],na.rm=TRUE))

Results_B <- data.frame(Statistic, AUCt_B, AUCinf_B, Cmax_B, Tmax_B, Lambda_B, Thalf_B, AUCRatio_B)

#Calculate Descriptive Statistics for Treatments A and B Combined:
AUCt_all <- c(length(na.omit(AUCt)),mean(AUCt,na.rm=TRUE),median(AUCt,na.rm=TRUE),min(AUCt,na.rm=TRUE),max(AUCt,na.rm=TRUE),sd(AUCt,na.rm=TRUE),100*sd(AUCt,na.rm=TRUE)/mean(AUCt,na.rm=TRUE))
AUCinf_all <- c(length(na.omit(AUCinf)),mean(AUCinf,na.rm=TRUE),median(AUCinf,na.rm=TRUE),min(AUCinf,na.rm=TRUE),max(AUCinf,na.rm=TRUE),sd(AUCinf,na.rm=TRUE),100*sd(AUCinf,na.rm=TRUE)/mean(AUCinf,na.rm=TRUE))
Cmax_all <- c(length(na.omit(Cmax)),mean(Cmax,na.rm=TRUE),median(Cmax,na.rm=TRUE),min(Cmax,na.rm=TRUE),max(Cmax,na.rm=TRUE),sd(Cmax,na.rm=TRUE),100*sd(Cmax,na.rm=TRUE)/mean(Cmax,na.rm=TRUE))
Tmax_all <- c(length(na.omit(Tmax)),mean(Tmax,na.rm=TRUE),median(Tmax,na.rm=TRUE),min(Tmax,na.rm=TRUE),max(Tmax,na.rm=TRUE),sd(Tmax,na.rm=TRUE),100*sd(Tmax,na.rm=TRUE)/mean(Tmax,na.rm=TRUE))
Lambda_all <- c(length(na.omit(Lambda)),mean(Lambda,na.rm=TRUE),median(Lambda,na.rm=TRUE),min(Lambda,na.rm=TRUE),max(Lambda,na.rm=TRUE),sd(Lambda,na.rm=TRUE),100*sd(Lambda,na.rm=TRUE)/mean(Lambda,na.rm=TRUE))
Thalf_all <- c(length(na.omit(T1.2)),mean(T1.2,na.rm=TRUE),median(T1.2,na.rm=TRUE),min(T1.2,na.rm=TRUE),max(T1.2,na.rm=TRUE),sd(T1.2,na.rm=TRUE),100*sd(T1.2,na.rm=TRUE)/mean(T1.2,na.rm=TRUE))
AUCRatio_all <- c(length(na.omit(AUCRatio)),mean(AUCRatio,na.rm=TRUE),median(AUCRatio,na.rm=TRUE),min(AUCRatio,na.rm=TRUE),max(AUCRatio,na.rm=TRUE),sd(AUCRatio,na.rm=TRUE),100*sd(AUCRatio,na.rm=TRUE)/mean(AUCRatio,na.rm=TRUE))

Results_all <- data.frame(Statistic, AUCt_all, AUCinf_all, Cmax_all, Tmax_all, Lambda_all, Thalf_all, AUCRatio_all)

#Report Descriptive Statistics for Treatment A:
Results_A
#Report Descriptive Statistics for Treatment B:
Results_B
#Report Descriptive Statistics for Treatments Combined:
Results_all

#Report ANOVA Results:
ANOVA_lnResults
ANOVA_Results