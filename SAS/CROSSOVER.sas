/* 2-Way Crossover Bioequivalence Program */

%LET PROJECT=TPD Part A Example;        /* project number */
%LET DRUG=Drug;           /* drug name */
%LET FLUID=Plasma;	      /* Plasma or Serum */
%LET PRFX=ng;             /* the prefix for "/mL" and ".h/mL". Usually mcg, ng, or pg. */
%LET TIMEFORMAT=8.2;	  /* format for Tmax, t1/2 */
%LET CONCFORMAT=9.2;	  /* Same as LLOQ. Format for Cmax, AUCt, etc. */
%LET RATIOFORMAT=7.4;	  /* format for AUCratio, lambda */

/* Output the results to a file. In SAS Studio, you need to create this file first
on the server, right-click the file, find out the folder, then copy the correct path and
filename below. */
FILENAME ZZ "/home/u63317948/output.txt";
PROC PRINTTO NEW PRINT=ZZ;
RUN;

OPTION NODATE PS=43 LS=120 NONUMBER;
PROC FORMAT;
  value $ParamFmt
			  'AUCt'="AUCt (&PRFX..h/mL)"
              'AUCinf'="AUCinf (&PRFX..h/mL)"
			  'Cmax'="Cmax (&PRFX/mL)"
			  'Tmax'="Tmax (h)"
			  'lambda'="Lambda (1/h)"
			  'halflife'="t1/2 (h)"
			  'AUCratio'="AUCt/AUCinf"
			  'lnauct'="AUCt (&PRFX..h/mL)"
			  'lnaucinf'="AUCinf (&PRFX..h/mL)"
			  'lncmax'="Cmax (&PRFX/mL)";
	VALUE PowerFmt 99.95-HIGH='  >99.9' LOW-99.9499999=[7.1];
RUN;
*Input the data.;
data avgBE2x2;
infile datalines DSD delimiter='09'x; /* tab delimited data */
input subject sequence $ period treatment $ AUCt AUCinf Cmax Tmax lambda halflife AUCratio;
lnauct=log(AUCt);
lnaucinf=log(AUCinf);
lncmax=log(Cmax);
datalines;
1	AB	1	A	364.74595	414.6804961	122.2	1.5	0.300192976	2.309005328	0.879583085
1	AB	2	B	375.426	422.7889664	126.2	1.5	0.266030635	2.605516397	0.887974923
2	BA	1	B	595.03875	612.1020981	206.9	1.5	0.315881735	2.194324976	0.972123363
2	BA	2	A	404.9456	437.9053151	102	1.5	0.250002161	2.77256476	0.924733238
3	BA	1	B	471.1644	499.6966333	122.8	1.5	0.220452424	3.144203032	0.942900889
3	BA	2	A	702.8314	759.9302729	201.5	0.66	0.255521681	2.712674626	0.92486301
4	AB	1	A	233.2503	259.9083832	59.47	3	0.328605771	2.109357903	0.897432769
4	AB	2	B	190.39375	.	37.26	1	.	.	.
5	BA	1	B	257.45585	284.6547633	84.67	2	0.3114095	2.225838262	0.904449471
5	BA	2	A	247.41105	261.7412831	66.4	1	0.429162593	1.615115558	0.94525039
6	AB	1	A	178.19155	210.1120445	54.19	1.5	0.261587426	2.649772546	0.848078702
6	AB	2	B	175.37495	189.7021219	55.27	1.5	0.543722099	1.274818849	0.924475426
7	BA	1	B	381.824	395.2178852	218.7	1	0.40466227	1.712902911	0.966110124
7	BA	2	A	246.3878	265.0736023	100.9	1	0.365518155	1.896341319	0.929507117
8	AB	1	A	407.99185	452.108046	89.51	1.5	0.171138962	4.050200909	0.902421122
8	AB	2	B	360.8275	406.4666669	181.9	0.66	0.405353587	1.709981614	0.887717319
9	BA	1	B	218.47035	241.7506936	59.68	1.5	0.298535113	2.321827987	0.903701027
9	BA	2	A	315.4761	383.042151	154.8	1.5	0.293342584	2.362927231	0.823606747
10	AB	1	A	140.1254	.	56.88	1	.	.	.
10	AB	2	B	91.80655	103.122699	25.56	2	0.485147377	1.428735295	0.890265198
11	AB	1	A	165.365	200.2484275	23.15	4	0.148494583	4.66782807	0.825799244
11	AB	2	B	269.01575	322.881255	57.05	1.5	0.141092152	4.912726701	0.8331724
12	BA	1	B	105.5625	125.9648343	47.2	0.66	0.356331774	1.945229784	0.838031507
12	BA	2	A	87.9882	112.2665638	37.76	0.66	0.262785419	2.6376927	0.783743592
13	BA	1	B	290.142	310.4973322	70.88	1.5	0.402842848	1.720639161	0.934442811
13	BA	2	A	182.77325	214.6103051	43.3	1	0.241228342	2.873406896	0.851651788
14	AB	1	A	122.47815	149.0976283	68.25	0.66	0.478596908	1.448290135	0.821462765
14	AB	2	B	230.48885	262.7033925	97.46	1.5	0.38926519	1.780655447	0.877372948
15	BA	1	B	143.5487	157.378009	88.38	1.5	0.461339029	1.502468112	0.912126802
15	BA	2	A	67.9815	.	27.54	1.5	.	.	.
16	AB	1	A	274.5791	296.1500222	60.43	2	0.254509286	2.723465189	0.927162179
16	AB	2	B	344.4828	370.9157081	98.82	2	0.263308145	2.632456283	0.928736078
;
run;
* Sort the data;
PROC SORT DATA=avgBE2x2;
  BY TREATMENT SUBJECT;
  RUN;
* Macro to report the data as inputted.;
%macro Repinput(data);
  options formdlim='';
  proc print label data=&data noobs;
  var subject period AUCt AUCinf Cmax Tmax lambda halflife AUCratio;
  BY Treatment;
  FORMAT AUCt AUCinf Cmax &CONCFORMAT;
  FORMAT Tmax halflife &TIMEFORMAT;
  FORMAT lambda AUCratio &RATIOFORMAT;
  title1 "&PROJECT - Individual &DRUG &FLUID Pharmacokinetic Parameters";
  title2 "(&SYSDATE9 - &SYSTIME)";
  label subject="Subject"
	  period="Period"
	  treatment="Treatment"
	  AUCt="AUCt (&PRFX..h/mL)"
      AUCinf="AUCinf (&PRFX..h/mL)"
	  Cmax="Cmax (&PRFX/mL)"
	  Tmax="Tmax (h)"
	  lambda="Lambda (1/h)"
	  halflife="t1/2 (h)"
	  AUCratio="AUCt/AUCinf";
run;
%mend Repinput;
%macro Repmeans(data);
/* To get rid of the parameter variable names in the output of PROC MEANS: */
ods noptitle; 
proc template;
   edit Base.Summary;
      column class nobs id type ways /*(varname)*/ (label) (min) (max) (range) (n)
        (nmiss) (sumwgt) (sum) (mean) (uss) (css) (var) (stddev) (cv) (stderr)
        (t) (probt) (lclm) (uclm) (skew) (kurt) (median) (mode) (q1)
        (q3) (qrange) (p1) (p5) (p10) (p25) (p50) (p75) (p90) (p95) (p99);
   end;
run;
options formdlim='';
title1 "&PROJECT - Descriptive Statistics for Individual &DRUG &FLUID Pharmacokinetic Parameters";
title2 "(&SYSDATE9 - &SYSTIME)";
PROC MEANS N MIN MAX MEDIAN MEAN STD CV MAXDEC=4 data=&data;
  VAR AUCt AUCinf Cmax Tmax lambda halflife AUCratio lnauct lnaucinf lncmax;
  BY treatment;
  label treatment="Treatment"
	  AUCt="AUCt (&PRFX..h/mL)"
      AUCinf="AUCinf (&PRFX..h/mL)"
	  Cmax="Cmax (&PRFX/mL)"
	  Tmax="Tmax (h)"
	  lambda="Lambda (1/h)"
	  halflife="t1/2 (h)"
	  AUCratio="AUCt/AUCinf"
      lnauct="ln(AUCt)"
	  lnaucinf="ln(AUCinf)"
	  lncmax="ln(Cmax)";*/
RUN;
%mend Repmeans;

*Macro to perform the ANOVA, calculate the 90%CI, and ISVs for ln-transformed parameters;
%macro lnci(data,var);
options formdlim='';
proc mixed data=&data;
  title1 "&PROJECT - ANOVA for ln-Transformed &DRUG &FLUID Pharmacokinetic Parameters";
  title2 "(&SYSDATE9 - &SYSTIME)";
  class subject sequence period treatment;
  model &var = sequence period treatment / ddfm=kenwardroger;
  random subject(sequence);
  lsmeans treatment / pdiff cl alpha=0.1;
  estimate "A vs B"  treatment 1 -1;
  ods output estimates=se_&var;
  make "CovParms" out=cov_&var;
  ods output "Least Squares Means"=LSM_&var;
run;
data LSM_&var; *Add the parameter name to the output file;
  set LSM_&var;
  Param="&var";
  length Param $8;
run;
* Now calculate the 90% confidence intervals;
options formdlim='';
data CI_&var;
  set se_&var;
  if _N_=1;
  Param="&var";
  length Param $8;
  Ratio=100*exp(ESTimate);
  LCI=100*exp(ESTimate-tinv(0.95,DF)*StdErr);
  UCI=100*exp(ESTimate+tinv(0.95,DF)*StdErr);
run;
* Now calculate the ISVs;
options formdlim='';
data ISV_&var;
   set cov_&var;
   Param="&var";
   length Param $8;
   IF covparm = "subject(sequence)" THEN CVtype="InterCV";
   else CVtype="IntraCV";
   ISV=100*sqrt(exp(estimate)-1);
 run;
%mend lnci;
* The main program is here;
%Repinput(avgBE2x2);
%Repmeans(avgBE2x2);
%lnci(avgBE2x2,lnauct);
%lnci(avgBE2x2,lnaucinf);
%lnci(avgBE2x2,lncmax);
*Merge data sets for printing in one appendix.;
DATA ALLCI;
  LENGTH Param $ 8.;
  SET CI_lnauct CI_lnaucinf CI_lncmax;
 RUN;
DATA ISV;
  LENGTH Param $ 8.;
  SET ISV_lnauct ISV_lnaucinf ISV_lncmax;
RUN;
DATA LSM_A;
  LENGTH Param $ 8.;
  SET LSM_lnauct LSM_lnaucinf LSM_lncmax;
  IF treatment="A";
  LSM_A = Estimate;
  GEO_A = exp(Estimate);
RUN;
DATA LSM_B;
  LENGTH Param $ 8.;
  SET LSM_lnauct LSM_lnaucinf LSM_lncmax;
  IF treatment="B";
  LSM_B = Estimate;
  GEO_B = exp(Estimate);
RUN;
* Merge data sets to print out one final appendix;
PROC SORT data=ISV;
by Param;
run;
PROC TRANSPOSE data=ISV out=ISVtransp let;
 BY Param;
 ID CVtype;
 VAR ISV;
run;
PROC SORT data=LSM_A;
by Param;
run;
PROC SORT data=LSM_B;
by Param;
run;
DATA LSM_both;
 MERGE LSM_A(KEEP=Param LSM_A GEO_A) LSM_B(KEEP=Param LSM_B GEO_B);
 BY Param;
RUN;
PROC SORT data=ISVtransp;
 BY Param;
run;
PROC SORT data=ALLCI;
 BY Param;
run;
DATA ALLCIandISV;
 MERGE ALLCI(KEEP=Param Probt Ratio LCI UCI StdErr DF) ISVtransp(KEEP=Param InterCV IntraCV);
 BY Param;
RUN;
DATA Summary;
 MERGE LSM_both ALLCIandISV(KEEP=Param Ratio Probt IntraCV InterCV LCI UCI StdErr DF);
 BY Param;
 *Now work out the power.;
 sigw=SQRT(LOG((IntraCV/100)**2+1));
 t1=tinv(0.95,DF);
 t2=-1*t1;
 tau2=(SQRT(DF+2)*(LOG(RATIO/100)-LOG(1.25)))/SQRT(2*sigw**2);
 tau1=(SQRT(DF+2)*(LOG(RATIO/100)-LOG(0.80)))/SQRT(2*sigw**2);
 probt2=probt(t2,DF,tau2);
 probt1=probt(t1,DF,tau1);
 POWER=100*(probt2-probt1);
RUN;
options formdlim='';
proc print label data=Summary noobs;
  var LSM_A LSM_B GEO_A GEO_B Probt Ratio LCI UCI IntraCV InterCV StdErr POWER;
  id Param;
  FORMAT Param $ParamFmt.;
  FORMAT POWER PowerFmt.;
  FORMAT LSM_A LSM_B GEO_A GEO_B StdErr &CONCFORMAT;
  FORMAT Ratio LCI UCI IntraCV InterCV 6.2;
  title1 "&PROJECT - Summary of Statistical Analysis of Ln-Transformed &DRUG &FLUID Data";
  title2 "(&SYSDATE9 - &SYSTIME)";
  label Param="Parameter"
	  LSM_A="Test Least Squares Mean"
      LSM_B="Reference Least Squares Mean"
	  GEO_A="Test Geometric Mean"
	  GEO_B="Reference Geometric Mean"
	  Ratio="Ratio A/B"
	  Probt="P-value of ANOVA"
	  IntraCV="Intra- Subject CV(%)"
	  InterCV="Inter- Subject CV(%)"
	  LCI="90% CI Lower Limit"
      UCI="90% CI Upper Limit"
      StdErr="Standard Error"
	  POWER="Power (%)";
run;
