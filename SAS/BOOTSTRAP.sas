/* 
Bootstrapping Program

This program was adapted from the programs accompanying the book:

Patterson S. and Jones B., Bioequivalence and Statistics in Clinical Pharmacology.
Chapman & Hall/CRC, Boca Raton FL, USA; 2006

with minor modifications. The original program is available for download here:
http://www.crcpress.com/e_products/downloads/download.asp?cat_no=C5300
*/

/* Output the results to a file. In SAS Studio, you need to create this file first
on the server, right-click the file, find out the folder, then copy the correct path and
filename below. */
FILENAME ZZ "/home/u63317948/output.txt";
PROC PRINTTO NEW PRINT=ZZ;
RUN;

data bootexample;
infile datalines DSD delimiter='09'x; /* tab delimited data */
input Subject Sequence$ Period formula$ AUC AUCinf Cmax;
lnauc=log(AUC);
lnaucinf=log(AUCinf);
lncmax=log(Cmax);
datalines;
1	TR	1	T	364.74595	414.6804961	122.2
1	TR	2	R	375.426	422.7889664	126.2
2	RT	1	R	595.03875	613.6228186	206.9
2	RT	2	T	404.9456	437.9053151	102
3	RT	1	R	471.1644	492.1550519	122.8
3	RT	2	T	702.8314	759.9302729	201.5
4	TR	1	T	233.2503	.	59.47
4	TR	2	R	190.39375	.	37.26
5	RT	1	R	257.45585	284.6547633	84.67
5	RT	2	T	247.41105	263.1714446	66.4
6	TR	1	T	178.19155	210.1120445	54.19
6	TR	2	R	175.37495	191.5725385	55.27
7	RT	1	R	381.824	395.2178852	218.7
7	RT	2	T	246.3878	265.0736023	100.9
8	TR	1	T	407.99185	437.4545188	89.51
8	TR	2	R	360.8275	406.4666669	181.9
9	RT	1	R	218.47035	237.8606236	59.68
9	RT	2	T	315.4761	383.042151	154.8
10	TR	1	T	140.1254	244.1144584	56.88
10	TR	2	R	91.80655	.	25.56
11	TR	1	T	165.365	200.2484275	23.15
11	TR	2	R	269.01575	322.881255	57.05
12	RT	1	R	105.5625	125.9648343	47.2
12	RT	2	T	87.9882	112.2665638	37.76
13	RT	1	R	290.142	310.4973322	70.88
13	RT	2	T	182.77325	214.6103051	43.3
14	TR	1	T	122.47815	149.0976283	68.25
14	TR	2	R	230.48885	264.9035429	97.46
15	RT	1	R	143.5487	157.378009	88.38
15	RT	2	T	67.9815	.	27.54
16	TR	1	T	274.5791	296.1500222	60.43
16	TR	2	R	344.4828	370.9157081	98.82
;
run;

ods listing;run;

***********  BOOTSTRAP MACRO;
%macro bootstrp(indata=,        /* dataset to sample from */
                seq=,           /* Sequence bootstrapping */
                seed=33920,     /* random number seed */
                nrep=,          /* # of bootstrap repetitions */
                nsamp=0,        /* # in bootstrap sample */
                bootsamp=);     /* output dataset */
   data test;
      set &indata;if sequence=&seq;
      run;
   data _null_;
      set test nobs=count;
      call symput('count',left(put(count,18.)));
      if &nsamp<=0 then call symput('nsamp',left(put(count,18.)));
      if &seed=0 then do;
         seed=time();
         call symput('seed',left(put(seed,18.)));
      end;
   run;
   data &bootsamp(label=seed=&seed);
      retain seed &seed;
      drop ijlkrst r seed;
      %do i=1 %to &nrep;
         rep=&i;
         put "generating rep " rep "of &nrep";
         do ijlkrst=1 to &nsamp;
            call ranuni(seed,r);
            pointvar = ceil(r*&count);
            set test point=pointvar;
            output;
         end;
      %end;
      stop;
   run;
proc sort data=&bootsamp;by rep;run;
%mend bootstrp;

option nomprint nomlogic nosymbolgen;

data all;set work.bootexample;run;

data pkt(keep=sequence subject T_auc T_cmax);
set all;
if formula='T' and sequence='TR' and period=1 then do;T_auc=auc;T_cmax=cmax;output;end;
if formula='T' and sequence='RT' and period=2 then do;T_auc=auc;T_cmax=cmax;output;end;
run; proc sort;by subject sequence;run;

data pkr(keep=sequence subject R_auc R_cmax);
set all;
if formula='R' and sequence='RT' and period=1 then do;R_auc=auc;R_cmax=cmax;output;end;
if formula='R' and sequence='TR' and period=2 then do;R_auc=auc;R_cmax=cmax;output;end;
run; proc sort;by subject sequence;run;

data pk(keep=sequence subject T_auc T_cmax R_auc R_cmax);
merge pkt pkr;by subject sequence;run;

%bootstrp(indata=pk,seq='TR',seed=43135,nrep=1000,nsamp=50,bootsamp=tr);run;
%bootstrp(indata=pk,seq='RT',seed=44837,nrep=1000,nsamp=50,bootsamp=rt);run;

data for_mix;
set tr rt;
by rep;
run;
data for_mix;set for_mix;
subj_id+1;
run;

data pkT(keep=rep sequence subject subj_id formula auc cmax);
set for_mix;
formula='T';
auc=T_auc;cmax=T_cmax;
run;proc sort;by rep subj_id sequence;run;
data pkR(keep=rep sequence subject subj_id formula auc cmax);
set for_mix;
formula='R';
auc=R_auc;cmax=R_cmax;
run;proc sort;by rep subj_id sequence;run;

data for_mix(keep=rep sequence subject subj_id formula auc cmax);
merge pkt pkr;by rep subj_id sequence formula;run;
run;

data for_auc(keep=rep sequence subject subj_id period formula lnauc);
set for_mix;
lnauc=log(auc);
if sequence='TR' then do;
        if formula='T' then do; period=1; end;
        if formula='R' then do; period=2; end;
        output;end;
if sequence='RT' then do;
        if formula='R' then do; period=1; end;
        if formula='T' then do; period=2; end;
        output;end;
run;
proc sort;by rep sequence subj_id period formula;run;

data for_cmax(keep=rep sequence subject subj_id period formula lncmax);
set for_mix;
lncmax=log(cmax);
if sequence='TR' then do;
        if formula='T' then do; period=1; end;
        if formula='R' then do; period=2; end;
        output;end;
if sequence='RT' then do;
        if formula='R' then do; period=1; end;
        if formula='T' then do; period=2; end;
        output;end;
run;
proc sort;by rep sequence subj_id period formula;run;

ods listing close;run;

proc mixed scoring=50 maxiter=200 data=for_auc method=reml;
by rep;*where rep <=1000;
class sequence subj_id period formula;
model lnauc=sequence period formula/ddfm=KENWARDROGER;
random int/subject=subj_id(sequence);
lsmeans formula;
estimate 'T-R' formula -1 1/cl alpha=0.10;
ods output Estimates=test_auc; 
run;

proc mixed scoring=50 maxiter=200 data=for_cmax method=reml;
by rep;
class sequence subj_id period formula;
model lncmax=sequence period formula/ddfm=KENWARDROGER;
random int/subject=subj_id(sequence);
lsmeans formula;
estimate 'T-R' formula -1 1/cl alpha=0.10;
ods output Estimates=test_cma; 
run;

ods listing;run;

proc sort data=test_auc;by rep;run;

data i_auc(keep=rep i_auc);set test_auc;
i_auc='F';
if upper<(log(1.25)) and lower>(log(0.8)) then i_auc='S';run;

proc sort data=test_cma;by rep;run;
data i_cmax(keep=rep i_cmax);set test_cma;
i_cmax='F';
if upper<(log(1.25)) and lower>(log(0.8)) then i_cmax='S';run;

data index;
merge i_auc i_cmax;by rep;
index='F';
if i_auc='S' and i_cmax='S' then index='S';
run;

title 'Number of Successful Data Based BE Simulations';
proc freq data=index;
tables i_auc i_cmax index;
run;

quit;run;
