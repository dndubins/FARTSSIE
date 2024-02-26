/* 
The procedure for running the ANOVA in this program was adapted from the book:

Patterson S. and Jones B., Bioequivalence and Statistics in Clinical Pharmacology.
Chapman & Hall/CRC, Boca Raton FL, USA; 2006

with minor modifications.
*/

/* Output the results to a file. In SAS Studio, you need to specify the correct directory on the server. 
Right-click on the SAS program, click "Properties", find out the folder location, then copy the correct 
path and output filename below. The program will create the output file, or over-write it if it exits.*/
FILENAME _n&sysindex "/home/u63317948/output/crossover.txt";
PROC PRINTTO NEW PRINT=_n&sysindex;
RUN;

OPTION NODATE PS=43 LS=120 NONUMBER; 
/* Input the data: */
data ABEexample1;
infile datalines DSD delimiter='09'x; /* tab delimited data */
input Subject$ Sequence$ Period Treatment$ AUCt AUCinf Cmax;
logauct=log(AUCt);
logaucinf=log(AUCinf);
logcmax=log(Cmax);
datalines;
1	AB	1	A	364.74595	414.6804961	122.2
2	BA	2	A	404.9456	437.9053151	102
3	BA	2	A	702.8314	759.9302729	201.5
4	AB	1	A	233.2503	259.9083832	59.47
5	BA	2	A	247.41105	261.7412831	66.4
6	AB	1	A	178.19155	210.1120445	54.19
7	BA	2	A	246.3878	265.0736023	100.9
8	AB	1	A	407.99185	452.108046	89.51
9	BA	2	A	315.4761	383.042151	154.8
10	AB	1	A	140.1254	.	56.88
11	AB	1	A	165.365	200.2484275	23.15
12	BA	2	A	87.9882	112.2665638	37.76
13	BA	2	A	182.77325	214.6103051	43.3
14	AB	1	A	122.47815	149.0976283	68.25
15	BA	2	A	67.9815	.	27.54
16	AB	1	A	274.5791	296.1500222	60.43
1	AB	2	B	375.426	422.7889664	126.2
2	BA	1	B	595.03875	612.1020981	206.9
3	BA	1	B	471.1644	499.6966333	122.8
4	AB	2	B	190.39375	.	37.26
5	BA	1	B	257.45585	284.6547633	84.67
6	AB	2	B	175.37495	189.7021219	55.27
7	BA	1	B	381.824	395.2178852	218.7
8	AB	2	B	360.8275	406.4666669	181.9
9	BA	1	B	218.47035	241.7506936	59.68
10	AB	2	B	91.80655	103.122699	25.56
11	AB	2	B	269.01575	322.881255	57.05
12	BA	1	B	105.5625	125.9648343	47.2
13	BA	1	B	290.142	310.4973322	70.88
14	AB	2	B	230.48885	262.7033925	97.46
15	BA	1	B	143.5487	157.378009	88.38
16	AB	2	B	344.4828	370.9157081	98.82
;
run;

/* Here's the ANOVA for AUCt: */
proc mixed data=ABEexample1;
class Sequence Subject Period Treatment;
model logauct=Sequence Period Treatment/ddfm=kenwardroger;
random Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for log(AUCt)' Treatment 1 -1;
run;

/* Here's the ANOVA for AUCinf: */
proc mixed data=ABEexample1;
class Sequence Subject Period Treatment;
model logaucinf=Sequence Period Treatment/ddfm=kenwardroger;
random Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for log(AUCinf)' Treatment 1 -1;
run;

/* Here's the ANOVA for Cmax: */
proc mixed data=ABEexample1;
class Sequence Subject Period Treatment;
model logcmax=Sequence Period Treatment/ddfm=kenwardroger;
random Subject(Sequence);
lsmeans Treatment/pdiff cl alpha=0.1;
estimate 'Average Bioequivalence for log(Cmax)' Treatment 1 -1;
run;
