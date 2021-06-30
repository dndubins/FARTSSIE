#1) Replace all missing values in the data set with "NA"
#2) Select the "Results" range from "Noncompartmental PK Analysis" (including headers), Control-C
#3) Now type in this line and hit Enter:
#For PC:
#ABEexample1 <- read.table(file="clipboard", header=TRUE)
#For MacOS:
#ABEexample1 <- read.table(pipe("pbpaste"), header=TRUE)

#Now copy and paste the following into R: (change variable as necessary)

library(outliers)

attach(ABEexample1)

ABE_trtA <- ABEexample1[Treatment == "A", ]
ABE_trtB <- ABEexample1[Treatment == "B", ]

grubbs.test(ABE_trtA$AUCt, type = 10, opposite = FALSE, two.sided = TRUE)
grubbs.test(ABE_trtB$AUCt, type = 10, opposite = FALSE, two.sided = TRUE)

indAUCratios <- ABE_trtA$AUCt / ABE_trtB$AUCt
grubbs.test(indAUCratios, type = 10, opposite = FALSE, two.sided = TRUE)

indAUCratios <- cbind(ABE_trtA$Subject,indAUCratios)
indAUCratios
