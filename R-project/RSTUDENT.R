#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers), Control-C
#3) Now type in this line and hit Enter:
#For PC:
#ABEexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)

factors       <- c("Subject", "Sequence", "Period", "Treatment")
ABEexample1[factors] <- lapply(ABEexample1[factors], factor) # factorize for glm()

#Now copy and paste the following into R: (change variable as necessary)
lm.lnauct <- glm(log(AUCt)~Treatment+Period+Sequence+Subject%in%Sequence, data=ABEexample1, na.action=na.exclude)

#to show residuals and r-studentized residuals:
rs.lnauct2 <- data.frame(Subject=ABEexample1$Subject,Period=ABEexample1$Period,Treatment=ABEexample1$Treatment,Residual=residuals(lm.lnauct),rStudent=rstudent(lm.lnauct))
rs.lnauct2

#Now copy and paste the following into R: (change variable as necessary)
lm.lnaucinf <- glm(log(AUCinf)~Treatment+Period+Sequence+Subject%in%Sequence, data=ABEexample1, na.action=na.exclude)

#to show residuals and r-studentized residuals:
rs.lnaucinf2 <- data.frame(Subject=ABEexample1$Subject,Period=ABEexample1$Period,Treatment=ABEexample1$Treatment,Residual=residuals(lm.lnaucinf),rStudent=rstudent(lm.lnaucinf))
rs.lnaucinf2

#Now copy and paste the following into R: (change variable as necessary)
lm.lncmax <- glm(log(Cmax)~Treatment+Period+Sequence+Subject%in%Sequence, data=ABEexample1, na.action=na.exclude)

#to show residuals and r-studentized residuals:
rs.lncmax2 <- data.frame(Subject=ABEexample1$Subject,Period=ABEexample1$Period,Treatment=ABEexample1$Treatment,Residual=residuals(lm.lncmax),rStudent=rstudent(lm.lncmax))
rs.lncmax2
