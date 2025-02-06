* This program calculates ICER given
* 50% censoring, Unif(0,1) delay distribution,
* KM estimated since a
* Written by P. Khudyakov
* Last update: April 20, 2016
***************************************************;
options nocenter  ps=500  ls=80;

%let numhzr=02;
%let numbhz=12;
libname yll "/udd/stlxu/ICER_mycode/DST/bhz&numbhz/hr&numhzr";

%include "/udd/stlxu/ICER_mycode/DST/bhz&numbhz/hr&numhzr/calc10.sas";

data sim; infile "/udd/stlxu/ICER_mycode/DST/simdata/bhz&numbhz._sim_&numhzr..txt";

input Z1 Z2 T D S O;
group=1-S;
*TS=0;
flag_TS=0;
run;

%macro cycle(sim=, Z1=, Z2=, group=, O=, n=, ii=, delay=, outdata=, avg=) ;

libname yll "/udd/stlxu/ICER_mycode/DST/bhz&numbhz/hr&numhzr";

%do i=1 %to &n;

%calc(sim=&sim, Z1=&Z1, Z2=&Z2, group=&group, O=&O, ii=&i, outdata=&outdata, delay=&delay, avg=&avg);

%end;
%mend;

************************************;
* With delay & fixed for all & 0.01;
************************************;
data delay;
   call streaminit(1234567);
   do i=1 to 10000000;
      TS=rand('UNIFORM');
      output;
   end;
run;
proc means; var TS; run;

data randit;
   call streaminit(7654321);
   do i=1 to 10000000;
      R=rand('UNIFORM');
      output;
   end;
run;
proc means; var R; run;

data sim01; merge sim delay randit; run;
proc means; var TS R; run;

data sim01; set sim01;
   if S=0 then TS=0;
        if R<0.1 then do;
   if TS GE T then delete;
        end;
        else TS=0;
run;
proc means; var TS R; run;

data sim02; set sim01; where TS>0; S=0; group=1; T=TS; TS=0; D=0; flag_TS=1; run;

data all; set sim01; * sim02; 
if T>10 then do; T=10; D=0; end; 
run;
proc means ; var TS T; where S=0; run;
proc means ; var TS T; where S=1; run;

*endsas;
options NOSYNTAXCHECK;
%cycle(sim=all, Z1=Z1, Z2=Z2, group=group, O=O,n=1000, outdata=delay_2_10, delay=1, avg=0.5) ;




