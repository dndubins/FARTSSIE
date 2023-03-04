/* r-Studentized Residuals Program 

The procedure for estimating studentized residuals was adapted from:

http://www.math.wustl.edu/~sawyer/s475f05/onetworeg.sas

with minor modifications. 

*/

/* Output the results to a file. In SAS Studio, you need to create this file first
on the server, right-click the file, find out the folder, then copy the correct path and
filename below. */
FILENAME ZZ "/home/u63317948/output.txt";
PROC PRINTTO NEW PRINT=ZZ;
RUN;

*Input the data.;
data ABEexample1;
infile datalines DSD delimiter='09'x; /* tab delimited data */
input subject sequence $ period treatment $ AUCt AUCinf Cmax;
lnauct=log(AUCt);
lnaucinf=log(AUCinf);
lncmax=log(Cmax);
datalines;
1	AB	1	A	364.74595	414.6804961	122.2
1	AB	2	B	375.426	422.7889664	126.2
2	BA	1	B	595.03875	613.6228186	206.9
2	BA	2	A	404.9456	437.9053151	102
3	BA	1	B	471.1644	492.1550519	122.8
3	BA	2	A	702.8314	759.9302729	201.5
4	AB	1	A	233.2503	.	59.47
4	AB	2	B	190.39375	.	37.26
5	BA	1	B	257.45585	284.6547633	84.67
5	BA	2	A	247.41105	263.1714446	66.4
6	AB	1	A	178.19155	210.1120445	54.19
6	AB	2	B	175.37495	191.5725385	55.27
7	BA	1	B	381.824	395.2178852	218.7
7	BA	2	A	246.3878	265.0736023	100.9
8	AB	1	A	407.99185	437.4545188	89.51
8	AB	2	B	360.8275	406.4666669	181.9
9	BA	1	B	218.47035	237.8606236	59.68
9	BA	2	A	315.4761	383.042151	154.8
10	AB	1	A	140.1254	.	56.88
10	AB	2	B	91.80655	.	25.56
11	AB	1	A	165.365	200.2484275	23.15
11	AB	2	B	269.01575	322.881255	57.05
12	BA	1	B	105.5625	125.9648343	47.2
12	BA	2	A	87.9882	112.2665638	37.76
13	BA	1	B	290.142	310.4973322	70.88
13	BA	2	A	182.77325	214.6103051	43.3
14	AB	1	A	122.47815	149.0976283	68.25
14	AB	2	B	230.48885	264.9035429	97.46
15	BA	1	B	143.5487	157.378009	88.38
15	BA	2	A	67.9815	.	27.54
16	AB	1	A	274.5791	296.1500222	60.43
16	AB	2	B	344.4828	370.9157081	98.82
;
run;

/* Here's the ANOVA for ln-AUCt: */
proc GLM data=ABEexample1;   *Note that GLM is used here;
class Sequence Subject Period Treatment;
model lnauct=Sequence Subject(Sequence) Period Treatment / SS3;
output out=Outstats_AUCt r=resid rstudent=rstudent cookd=cookd;
Test H=Sequence E=Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for ln(AUCt)' Treatment 1 -1;
contrast 'A vs. B' Treatment 1 -1;
run;

/* Sort the data by treatment (required for PROC PRINT) */
PROC SORT DATA=Outstats_AUCt;
by Treatment Subject;
run;

/* Print the outlier statistics for AUCt */
options formdlim='';
proc print label data=Outstats_AUCt noobs;
var Subject Sequence Period AUCt lnauct resid rstudent cookd;
BY Treatment;
title1 "Outlier Analysis for ln-AUCt";
label subject="Subject"
	period="Period"
	treatment="Treatment"
	AUCt="AUCt"
    lnauct="ln(AUCt)";
run;


/* Here's the ANOVA for ln-AUCinf: */
proc GLM data=ABEexample1;   *Note that GLM is used here;
class Sequence Subject Period Treatment;
model lnaucinf=Sequence Subject(Sequence) Period Treatment / SS3;
output out=Outstats_AUCinf r=resid rstudent=rstudent cookd=cookd;
Test H=Sequence E=Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for ln(AUCinf)' Treatment 1 -1;
contrast 'A vs. B' Treatment 1 -1;
run;

/* Sort the data by treatment (required for PROC PRINT) */
PROC SORT DATA=Outstats_AUCinf;
by Treatment Subject;
run;

/* Print the outlier statistics for AUCinf */
options formdlim='';
proc print label data=Outstats_AUCinf noobs;
var Subject Sequence Period AUCinf lnaucinf resid rstudent cookd;
BY Treatment;
title1 "Outlier Analysis for ln-AUCinf";
label subject="Subject"
	period="Period"
	treatment="Treatment"
	AUCinf="AUCinf"
    lnaucinf="ln(AUCinf)";
run; 

/* Here's the ANOVA for ln-Cmax: */
proc GLM data=ABEexample1;   *Note that GLM is used here;
class Sequence Subject Period Treatment;
model lncmax=Sequence Subject(Sequence) Period Treatment / SS3;
output out=Outstats_Cmax r=resid rstudent=rstudent cookd=cookd;
Test H=Sequence E=Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for ln(Cmax)' Treatment 1 -1;
contrast 'A vs. B' Treatment 1 -1;
run;

/* Sort the data by treatment (required for PROC PRINT) */
PROC SORT DATA=Outstats_Cmax;
by Treatment Subject;
run;

/* Print the outlier statistics for Cmax */
options formdlim='';
proc print label data=Outstats_Cmax noobs;
var Subject Sequence Period Cmax lncmax resid rstudent cookd;
BY Treatment;
title1 "Outlier Analysis for ln-Cmax";
label subject="Subject"
	period="Period"
	treatment="Treatment"
	Cmax="Cmax"
    lncmax="ln(Cmax)";
run; 

