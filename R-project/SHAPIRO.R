#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers), Control-C
#3) Now type in this line and hit Enter:
#For PC:
#ABEexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)
#Now copy and paste the following into R: (change variable as necessary)

attach(ABEexample1)

ABE_trtA <- ABEexample1[Treatment == "A", ]
ABE_trtB <- ABEexample1[Treatment == "B", ]

shapiro.test(log(ABE_trtA$AUCt))
shapiro.test(log(ABE_trtB$AUCt))

shapiro.test(log(ABE_trtA$AUCinf))
shapiro.test(log(ABE_trtB$AUCinf))

shapiro.test(log(ABE_trtA$Cmax))
shapiro.test(log(ABE_trtB$Cmax))



