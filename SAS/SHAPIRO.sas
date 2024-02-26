/* Shapiro-Wilke Log Normality Program */

/* Output the results to a file. In SAS Studio, you need to specify the correct directory on the server. 
Right-click on the SAS program, click "Properties", find out the folder location, then copy the correct 
path and output filename below. The program will create the output file, or over-write it if it exits.*/
FILENAME _n&sysindex "/home/u63317948/output/shapiro.txt";
PROC PRINTTO NEW PRINT=_n&sysindex;
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

/* Sort the data by treatment (required for PROC UNIVARIATE) */
PROC SORT DATA=ABEexample1;
by Treatment Subject;
run;

/* Now we test the normality of the distribution for ln-AUCt. If the Shapiro-Wilk p-value is less than 0.05,
we reject the null hypothesis that the sample came from a normally-distributed population. */
title1 'PROC UNIVARIATE: Checking ln-AUCt';
title2 'for fit to a normal distribution';
proc univariate data=ABEexample1 normaltest;
by Treatment;
var lnauct;
run;

title1 'PROC UNIVARIATE: Checking ln-AUCinf';
title2 'for fit to a normal distribution';
proc univariate data=ABEexample1 normaltest;
by Treatment;
var lnaucinf;
run;

title1 'PROC UNIVARIATE: Checking ln-Cmax';
title2 'for fit to a normal distribution';
proc univariate data=ABEexample1 normaltest;
by Treatment;
var lncmax;
run;
