*options nocenter ps=500 ls=80; %macro hr(sim=, Z1=, Z2=, group=, O=, ii=,outhazard=) ;

libname yll "/udd/stlxu/ICER_mycode/DST/bhz&numbhz/hr&numhzr";

****************************************************;
* Find Cox PH fit;
****************************************************;
proc phreg data=&sim ;
  model (TS,T)*d(0)= &group  &Z1  / rl;
*  strata &group;
where &O=&ii;
ods output ParameterEstimates =hazard;
run;

data yll.&outhazard; set hazard; run;

%mend hr;


%macro calc(sim=, Z1=, Z2=, group=, O=, ii=,outdata=, delay=, avg=) ; *,TSf=, TSp=,fixed=,TSd=, avg=,);

**************************************************;
* Macro parameters;
**************************************************;
* delay - and indicator of delay;
* outdata - name of file with the output data;
* TSf - percent that delay takes from the total time;
* TSp - percent of subjects in the intervention group that have a delay;
* fixed - an indicator that delay is fixed for all subjects;
* TSd - value of the fixed delay;
**************************************************;
  
libname yll "/udd/stlxu/ICER_mycode/DST/bhz&numbhz/hr&numhzr";

data ETS; set yll.&outdata; run;
****************************************************;
* Find beta values for each x;
****************************************************;
proc phreg data=&sim covout outest=cov noprint;
  model (TS,T)*d(0)= &Z1 ;     
  strata &group;
where &O=&ii;
 output out=xb xbeta=xb      ;
run;

proc sort data=xb;  by T; run;

data xb; set xb;
ebz=exp(xb);
run;

proc sort data=xb; by &group ; run;

********************************;
* tp - times of events;
********************************;
data tp; set xb ;
  if d=1;
keep &group t;
run;
proc sort data=tp ;by &group t; run;

*********************************;
* np - number of events per group;
*********************************;
proc  means data=tp noprint ; by &group;
output out=np n=n ;
run;

data NP;set NP; keep n; run;

proc means data=xb  noprint; by &group;
output out=ng n=n ;
run;
data NG;set NG; keep n; run;

DATA ts; set xb; keep ts; run;
DATA x; set xb ;    keep t ebz; run;
data x; merge x ts; run;
proc sort data=xb; by &group t; run;

title 'X data';
proc print data=x (obs=10) ; run;
title ' ';

DATA z; set xb ;    keep &Z1 ; *&Z2 ; run;

data cov; set cov ; if _type_='COV  ' ; keep &Z1; * &Z2; *z1 -z8; run;

********************************************************************;
proc IML;
********************************************************************;
use tp ; read all var { t } into tp ;
use NP ; read all into np ;
use NG ; read all into ng ;
use  X ; read all into x ;
use  Z ; read all into z ;
use cov; read all into sig;
ND=10; 
m=1:10; 
m=m/10-0.05;

use ETS; read all into ets;

* &&&&&&&&&&&&&&&&&&  ;
          tu=10;               *    ++++++++++++++++++++++++++++;
* &&&&&&&&&&&&&&&&&&  ;
* Number of covariates;
nz=ncol(z) ;
* Number of groups;
nrg=nrow(np) ;
* Maximum of the number of events in groups;
nrt=max(np) ;
* Number of individuals;
nn=nrow(x);

s0=j(nrg,nrt,0) ;
ss0=j(nrg,nrt,0) ;
fh=j(nrg,nrt,0) ;
lh=j(nrg,nrt,0) ;
gm=j(nrg,nrt,0) ;
vh=j(nrg,nrt,0) ;
li=j(nrg,1,0) ;
h =j(nz ,nrt,0) ;

* Conditional measures;
**********************************************;
hC=j(nz*ND,nrt,0) ;
***********************************************;
m0=0.5; 
aaa=j(nrg,1,0);
a1=j(nrg,1,0);
a2=j(nrg,1,0);
ttt=j(nrg,1,0);
a12=j(nrg*ND,2,0);
**************************************;
* index for a in both groups;
**************************************;

do g=1 to nrg;
    do k=1 to np[1];
      if tp[k]<m0 then aaa[1]=k;
      do nn=1 to ND;
      if tp[k]<m[nn] then a12[nn,1]=k;
      end;
    end;
    do k=np[1]+1 to (np[1]+np[2]);
      if tp[k]<m0 then aaa[2]=k-np[1];
      do nn=1 to ND; 
      if tp[k]<m[nn] then a12[nn,2]=k-np[1];
      end;
    end;
end;  
**************************;
* Indexes for time vectors;
**************************;
p1=1; p2=0;
j1=1; j2=0;
  do   g=1 to nrow(ng) ;
p2=p2+np[g] ;
j2=j2+ng[g] ;
  k=1;
s1=j(nz,nrt,0) ;
SS=0;
do   p=p1 to p2    ;
  If (tp[p] <= tu) then li[g]= li[g] +1;

**************************************;
* find S0(t) and S1(t) for all t<=t[p];
**************************************;
  do    j=j1 to j2 ;
   if &delay=1 then do;
	if g=1 then
	s0[g,k]= s0[g,k] + x[j,2]*( x[j,1] >= tp[p] )*(x[j,3] <= tp[p] ) ;
        else s0[g,k]= s0[g,k] + x[j,2]*( x[j,1] >= tp[p] ) ;
   end;
   else s0[g,k]= s0[g,k] + x[j,2]*( x[j,1] >= tp[p] ) ;
   
         do  l=1 to nz;
   if &delay=1 then do;
        if g=1 then
         s1[l,k]= s1[l,k] + x[j,2]*z[j,l]*(x[j,1] >= tp[p] )*(x[j,3] <= tp[p] ) ;
	 else  s1[l,k]= s1[l,k] + x[j,2]*z[j,l]*(x[j,1] >= tp[p] ) ;
   end;
   else  s1[l,k]= s1[l,k] + x[j,2]*z[j,l]*(x[j,1] >= tp[p] ) ;
         end;
  end;
   SS=SS+1/s0[g,k];
   ss0[g,k]=exp(-SS);
   s0[g,k]=s0[g,k]/ng[g];
  do l=1 to nz;
    s1[l,k]=s1[l,k]/ng[g];
  end;
**********************************;
* Set values for the first event;
**********************************;
   IF k=1 then do;
        lh[g,k]= (1.0/s0[g,k])/ng[g] ;
        vh[g,k]=(1.0/s0[g,k]**2 )/ng[g];
   do l=1 to nz;
    h[l,k]= (s1[l,k]/s0[g,k]**2 )/ng[g];
   end;
*******************************;
*   if k>aaa[g] then do;
*          vhC[g,k]=(1.0/s0[g,k]**2 )/ng[g];
*        do l=1 to nz;
*          hC[l,k]=(s1[l,k]/s0[g,k]**2 )/ng[g];
*        end;
*    end;
********************************;
end;
else;  do;
        lh[g,k]=  lh[g,k-1] +(1.0/s0[g,k])/ng[g] ;
        vh[g,k]= vh[g,k-1] +(1.0/s0[g,k]**2 )/ng[g];
        do l=1 to nz;
          h[l,k]=  h[l,k-1] +(s1[l,k]/s0[g,k]**2 )/ng[g];
        end;
******************************************************************;  
 do d=1 to ND;
   if k>=a12[d,g] then do;
       do l=1 to nz;
          hC[(l+(d-1)),k]=hC[(l+(d-1)),k-1] +(s1[l,k]/s0[g,k]**2 )/ng[g];
       end;
   end;
  end;

*****************************************************************;
 end;
  k=k+1;
end;
 p1=p2+1;  j1=j2+1;

if g=1 then do; sa1=s1;  ha=h; haC=hC;  end;
if g=2 then do; sb1=s1;  hb=h; hbC=hC;  end;
end;
*****************************************************************;

* Conditional measures;
**********************************************;
fhC1=j(nrg*ND,nrt,0) ;
vhC1=j(nrg*ND,nrt,0) ;
gmC1=j(nrg*ND,nrt,0) ;

******************************************************************;
p1=1; p2=0;
j1=1; j2=0;

  do   g=1 to nrow(ng) ;

p2=p2+np[g] ;
j2=j2+ng[g] ;
  k=1;

chi=j(nz,nrt,0) ;
chiC=j(nz*ND,nrt,0) ;

do   p=p1 to p2    ;

*****************************************************************;
 do j = 1 to nrow(x) ;
        e1=exp(-lh[g,k]*x[j,2] ) ;
        fh[g,k]=fh[g,k] + e1;
        gm[g,k]=gm[g,k] +x[j,2]*e1 ;
    do l=1 to nz;
          chi[l,k]=chi[l,k] +z[j,l]*x[j,2]*e1/ng[+] ;
    end;
**********************************;
* By delay times;
**********************************;
  do d=1 to ND;
   if k>=a12[d,g] then do;
        e1C=exp(-(lh[g,k]-lh[g,a12[d,g]]+lh[2,a12[d,2]])*x[j,2]) ;
        fhC1[(g+(d-1)*nrg),k]=fhC1[(g+(d-1)*nrg),k] + e1C;
        gmC1[(g+(d-1)*nrg),k]=gmC1[(g+(d-1)*nrg),k] +x[j,2]*e1C ;
        vhC1[(g+(d-1)*nrg),k]= vhC1[(g+(d-1)*nrg),k-1] +(1.0/s0[g,k]**2 )/ng[g];
       do l=1 to nz;
	 chiC[(l+(d-1)),k]=chiC[(l+(d-1)),k] +z[j,l]*x[j,2]*e1C/ng[+] ;
       end;
   end;
  end;
**************************************;
* end of by delay loop;
*************************************;
 end;
   fh[g,k]=fh[g,k]/ ng[+]  ;
   gm[g,k]=gm[g,k]/ng[+] ;

**********************************;
* By delay times;
**********************************;
  do d=1 to ND;
   fhC1[(g+(d-1)*nrg),k]=fhC1[(g+(d-1)*nrg),k]/ ng[+];
   gmC1[(g+(d-1)*nrg),k]=gmC1[(g+(d-1)*nrg),k]/ng[+];
  end;
**************************************;
* end of by delay loop;
*************************************;

  k=k+1;
end;
 p1=p2+1;  j1=j2+1;
if g=1 then do; chia=chi;  chiaC=chiC;  end;
if g=2 then do; chib=chi;  chibC=chiC;  end;
end;

****************************************************************;
   ah=j(nrg,1,0);
   ah0=j(nrg,1,0);
   omga=j(nrg,1,0);
   omga00=j(nrg,1,0);
   omga01=j(nrg,1,0);

   Pa0_1=j(nrg,ND,0); 
   ah_1=j(nrg,ND,0);
   ah0_1=j(nrg,ND,0);
   ah00_1=j(nrg,ND,0);
   ah01_1=j(nrg,ND,0);
   omga_1=j(nrg,ND,0);
   omgaF_1=j(nrg,ND,0);
   omga00_1=j(nrg,ND,0);
   omga01_1=j(nrg,ND,0);

/*   Pa0_2=j(nrg,1,0);
   ah_2=j(nrg,1,0);
   ah0_2=j(nrg,1,0);
   ah00_2=j(nrg,1,0);
   ah01_2=j(nrg,1,0);
   omga_2=j(nrg,1,0);
   omgaF_2=j(nrg,1,0);
   omga00_2=j(nrg,1,0);
   omga01_2=j(nrg,1,0);
*/
****************************************************************;
   phia=j(nz,1,0);
   phib=j(nz,1,0);
   phiaF=j(nz,1,0);
   phibF=j(nz,1,0);
***********************************************;
   ww=j(nrg,ND,0);
   msum=sum(1/m);

	do kk=1 to ND;
	   ww[kk]=(1/m[kk])/msum;
	end;
***********************************************;

p1=2 ; p2=0;
  do g=1 to nrg; 

   phi_1=j(nz,ND,0);
   phi0_1=j(nz,ND,0);
   phiF_1=j(nz,ND,0);

*   phi_2=j(nz,1,0);
*   phi0_2=j(nz,1,0);
*   phiF_2=j(nz,1,0);
   chiC1=j(nz,nrt,0);


k=2;
    do  p=p1 to p2+li[g]+1 ;
 tpp= tu;
if k < li[g]+1 then tpp=tp[p] ;

* Scenario I;
*********************************************************************;
    do d=1 to ND;    
        if g=1 then do; h=ha; chi=chia; nz1=1+(d-1)*nz; nz2=d*nz; 
	hC1=haC[nz1:nz2,1:nrt]; chiC1=chiaC[nz1:nz2,1:nrt];
	end;
        if g=2 then do; h=hb; chi=chib; nz1=1+(d-1)*nz; nz2=d*nz; 
        hC1=hbC[nz1:nz2,1:nrt]; chiC1=chibC[nz1:nz2,1:nrt];
        end;

	if tp[p-1]<m[d] then do;
       * Pa0_1[g,d]=fh[g,k-1];
        ah0_1[g,d]=ww[d]*fhC1[(g+(d-1)*nrg),k-1]*(tpp-tp[p-1]); 
        ah00_1[g,d]=ah00_1[g,d]+ww[d]*fh[g,k-1]*(tpp-tp[p-1]);
        
	omga_1[g,d]=(ww[d])*vhC1[(g+(d-1)*nrg),k-1]*gmC1[(g+(d-1)*nrg),k-1]**2*(tpp-tp[p-1])**2;
***********;
        if g=2 then omga00_1[g,d]= omga00_1[g,d]+(ww[d])*vh[g,k-1]*gm[g,k-1]**2*(tpp-tp[p-1])**2;
        if g=1 then omga00_1[g,d]= omga00_1[g,d]+(ww[d])*vh[g,k-1]*gm[g,k-1]**2*(tpp-tp[p-1])**2;
***********;

        end;
   
	if tp[p-1]>=m[d] then do; 
	ah0_1[g,d]= ah0_1[g,d]+ww[d]*fhC1[(g+(d-1)*nrg),k-1]*(tpp-tp[p-1] ) ;
        ah01_1[g,d]= ah01_1[g,d]+ww[d]*fhC1[(g+(d-1)*nrg),k-1]*(tpp-tp[p-1] ) ;
	omga_1[g,d]= omga_1[g,d]+(ww[d])*vhC1[(g+(d-1)*nrg),k-1]*gmC1[(g+(d-1)*2),k-1]**2*(tpp-tp[p-1])**2;
***********;
	if g=1 then omga01_1[g,d]= omga01_1[g,d]+(ww[d])*vhC1[(g+(d-1)*nrg),k-1]*gmC1[(g+(d-1)*nrg),k-1]**2*(tpp-tp[p-1])**2;
	if g=2 then omga01_1[g,d]= omga01_1[g,d]+(ww[d])*vh[g,k-1]*gm[g,k-1]**2*(tpp-tp[p-1])**2;
***********;
	end;
	
	if tp[p-1]>m[d] then do;
           do l=1 to nz;
	phi_1[l,d]= phi_1[l,d] +(gmC1[(g+(d-1)*2),k-1]*(h[l,k-1])-(lh[g,k-1]-lh[g,a12[d,g]]+lh[2,a12[d,2]])*chiC1[l,k-1])*(tpp - tp[p-1]);
        phiF_1[l,d]=phiF_1[l,d]+(gmC1[(g+(d-1)*2),k-1]*(h[l,k-1])-(lh[g,k-1]-lh[g,a12[d,g]]+lh[2,a12[d,2]])*chiC1[l,k-1])*(tpp - tp[p-1]);
           end;
	end;
        else do;
           do l=1 to nz;
       phi0_1[l,d]= phi0_1[l,d] +(gm[g,k-1]*h[l,k-1]-lh[g,k-1]*chi[l,k-1])*(tpp - tp[p-1]);
           end;
        end;

    end;
 
  k=k+1;
  end;

p2=p2 +np[g] ;
p1=p2+2;
  do d=1 to ND; 
	if g=1 then do; 
            do l=1 to nz;
		phia[l]=phia[l]+ww[d]*phi_1[l,d];    *+ww[d]*phi_1[l,2]; 
		phiaF[l]=phiaF[l]+ww[d]*phiF_1[l,d]; *ww[d]+ww[d]*phiF_1[l,2]; 
            end;
        end;
	if g=2 then do; 
            do l=1 to nz;
		phib[l]= phib[l]+ww[d]*phi_1[l,d];                      *0.5+0.5*phi_1[l,2]; 
		phibF[l]=phibF[l]+ww[d]*phiF_1[l,d]+ww[d]*phi0_1[l,d];  *0.5+0.5*phiF_1[l,2]+0.5*phi0_1[l,1]+0.5*phi0_1[l,2]; 
                phiaF[l]=phiaF[l]+ww[d]*phi0_1[l,d];                    *0.5+0.5*phi0_1[l,2]; 
             end;
         end;
   end;
end;
************************************************************************************;              
 p1=2; p2=0; 
     do g=1 to nrg ;
 k=2;
       do p=p1 to (p2+li[g] )  ;
       kq=k+1;
            do q=(p+1) to (p2+li[g]+1) ;
            tpq= tu;
		if kq<li[g]+1 then tpq=tp[q] ;
 
	do d=1 to ND; 

        if g=1 then do;
                if (tp[p-1]>m[d]) then do;
                        if (tpq>m[d]) then do;
       omga_1[g,d]=omga_1[g,d]+(ww[d])*vhC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),kq-1]*(tp[p]-tp[p-1])*(tpq-tp[q-1])*2.0 ;
       omga01_1[g,d]=omga01_1[g,d]+(ww[d])*vhC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),kq-1]*(tp[p]-tp[p-1])*(tpq-tp[q-1])*2.0;
                        end;
                end;

	end;
	if g=2 then do;
                omga01_1[g,d]=omga01_1[g,d]+(ww[d])*vh[g,k -1]*gm[g,k-1]*gm[g,kq-1]*(tp[p]-tp[p-1])*(tpq-tp[q-1])*2.0 ;

		if (tp[p-1]>m[d]) then do;
			if (tpq>m[d]) then do;
                        omga_1[g,d]=omga_1[g,d]+(ww[d])*vhC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),k-1]*gmC1[(g+(d-1)*2),kq-1]*(tp[p]-tp[p-1])*(tpq -tp[q-1])*2.0 ;
	                end;
		end;
                else if (tp[p-1]<=m[d]) then do;
                   if (tpq<=m[d]) then do;
                   omga01_1[1,d]=omga01_1[1,d]+(ww[d])*vh[g,k-1]*gm[g,k-1]*gm[g,kq-1]*(tp[p]-tp[p-1])*(tpq-tp[q-1])*2.0 ;
                   end;
		end;
	end;
end;		
 		 kq=kq+1;
    	    end ;    
       k =k +1;     
       end;

     p2=p2 +np[g] ;
     p1=p2+2;
     end;

*******************************************************************************;     
* new; 
**********;
 p1=2; p2=0;
 k=2;
       do p=2 to li[1]  ;
       kq=np[1]+1;

	       do q=(np[1]+2) to (np[1]+li[2]+1) ;
       	       tpq= tu;
                 if kq< li[2]+1 then tpq=tp[q] ;
 	do d=1 to ND;
                  if (tp[p-1]>m[d]) then do;
                        if (tpq<=m[d]) then do;
                        omga01_1[1,d]=omga01_1[1,d]+(ww[d])*vh[2,kq-1]*gm[2,kq-1]*gmC1[(1+(d-1)*2),k-1]*(tpq -tp[q-1])*(tp[p]-tp[p-1])*2.0 ;
			end;
                  end;
            do dd=1 to ND;
              if d^=dd then do;
                  if (tp[p-1]>m[d]) then do;
                        if (tpq<=m[dd]) then do;
                       omga01_1[1,d]=omga01_1[1,d]+(ww[d]*ww[dd])*vh[2,kq-1]*gm[2,kq-1]*gmC1[(1+(d-1)*2),k-1]*(tpq -tp[q-1])*(tp[p]-tp[p-1]); 
                        end;
                  end;
                  if (tp[p-1]>m[dd]) then do;
                        if (tpq<=m[d]) then do;
                        omga01_1[1,dd]=omga01_1[1,dd]+(ww[d]*ww[dd])*vh[2,kq-1]*gm[2,kq-1]*gmC1[(1+(dd-1)*2),k-1]*(tpq -tp[q-1])*(tp[p]-tp[p-1]);
                        end;
                   end;
               end;
             end;
	end;
               kq=kq+1;
               end ;
       k =k +1;
       end;

*********************************************************************************;
/*	do i=1 to nrg;
        if Pa0_1[i,1]=0 then Pa0_1[i,1]=1;
        ah0_1[i,1]=ah0_1[i,1];
        if Pa0_1[i,2]=0 then Pa0_1[i,2]=1;
        ah0_1[i,2]=ah0_1[i,2];
        end;
*/
	do g=1 to nrg;
           do d=1 to ND;
	ah[g]=ah[g]+ah00_1[2,d]+ah01_1[g,d];
	ah0[g]=ah0[g]+ah0_1[g,d];

        omga[g]=omga[g]+omga_1[g,d];
	omga00[g]=omga00[g]+omga00_1[g,d]; 
	omga01[g]=omga01[g]+omga01_1[g,d];
	   end;
	end;

*omga[1]=omga_1[1,1]+omga_1[1,2];
*omga[2]=omga_1[2,1]+omga_1[2,2];

 v=(phia-phib)`*sig*(phia-phib)+omga[1]/(ng[1])+omga[2]/(ng[2]) ;

*ah[1]=(ah00_1[2,1]+ah00_1[2,2])+(ah01_1[1,1]+ah01_1[1,2]);
*ah[2]=(ah00_1[2,1]+ah00_1[2,2])+(ah01_1[2,1]+ah01_1[2,2]);

*ah0[1]=(ah0_1[1,1]+ah0_1[1,2]);
*ah0[2]=(ah0_1[2,1]+ah0_1[2,2]);

czi=(ah[1]-ah[2])/sqrt(v);
******************;
* cov(mu_0,mu_1);
******************;
*omga00[1]=omga00_1[1,1]+omga00_1[1,2];
*omga00[2]=omga00_1[2,1]+omga00_1[2,2];

*omga01[1]=omga01_1[1,1]+omga01_1[1,2];
*omga01[2]=omga01_1[2,1]+omga01_1[2,2];

va1=(omga00[2]+omga01[1])/(ng[1])+phiaF`*sig*phiaF;
vb1=(omga01[2])/(ng[2])+phibF`*sig*phibF;

v1=(omga00[2]+omga01[1])/(ng[1])+(omga00[2]+omga01[2])/(ng[2])+(phiaF - phibF)`*sig*(phiaF- phibF);

va=omga[1]/(ng[1])+phia`*sig*phia;
vb=omga[2]/(ng[2])+phib`*sig*phib;

cc=(va+vb-v)/2;
cc1=0; *(va1+vb1-v1)/2;
cc_check=0; *phia`*sig*phib;
****************************;
* Add costs;
****************************;
Ca=330; 
Cb=115;  

dTime=ah0[1]-ah0[2]; 
Costs=Ca*ah0[1]-Cb*ah0[2];
ICER=Costs/dTime;

nrgn=nrg+1;

grad=j(nrg,1,0);
grad[1]=(Cb-Ca)*ah0[2]/(ah0[1]-ah0[2])**2;
grad[2]=(Ca-Cb)*ah0[1]/(ah0[1]-ah0[2])**2;

SS=j(nrg,nrg,0);

SS[1,1]=va;
SS[2,2]=vb;
SS[1,2]=cc;
SS[2,1]=cc;

V_ICER=grad`*SS*grad;


*  ------------- ;
*iii=1;
/*ETS[&ii,1]=omga_1[1]; 
ETS[&ii,2]=omga_1[2];
ETS[&ii,3]=omga_2[1];  
ETS[&ii,4]=omga_2[2]; 
ETS[&ii,5]=omga00[2];   
ETS[&ii,6]=omga00_1[2];
ETS[&ii,7]=omga00_2[2];  
ETS[&ii,8]=phiaF`*sig*phiaF;
ETS[&ii,9]=phibF`*sig*phibF; 
ETS[&ii,10]=v; 
ETS[&ii,11]=phia`*sig*phia; 
ETS[&ii,12]=phib`*sig*phib;
ETS[&ii,13]=omga01_1[1];  
ETS[&ii,14]=omga01_1[2];  
ETS[&ii,15]=omga01_2[1];  
ETS[&ii,16]=omga01_2[2];  
ETS[&ii,17]=va;
ETS[&ii,18]=vb;
ETS[&ii,19]=va1;
ETS[&ii,20]=vb1;


ccc={omga11,omga12,omga21,omga22,omga0,omga01,omga02,Fa,Fb,v,pha,phb,omga_11,omga_12,omga_21,omga_22,va,vb,va1,vb1};
*/
ETS[&ii,1]=ah0[1];
ETS[&ii,2]=ah0[2];
ETS[&ii,3]=m0;
ETS[&ii,4]=ICER;
ETS[&ii,5]=V_ICER;
ETS[&ii,6]=va;
ETS[&ii,7]=vb;
ETS[&ii,8]=ah[1];
ETS[&ii,9]=ah[2];
ETS[&ii,10]=omga01_1[1,1];
ETS[&ii,11]=omga01_1[1,2];
ETS[&ii,12]=omga01_1[2,1];
ETS[&ii,13]=omga01_1[2,2];
ETS[&ii,14]=v;
ETS[&ii,15]=omga00[1];
ETS[&ii,16]=omga00[2];
ETS[&ii,17]=cc;
ETS[&ii,18]=cc1;
ETS[&ii,19]=va1;
ETS[&ii,20]=vb1;


ccc={ah01,ah02,m0,ICER,V_ICER,va,vb,ah1,ah2,omga01_1,omga01_2,omga02_1,omga02_2,v,omga00_1,omga00_2,cc,cc1,va1,vb1};

*ah001,ah002,ah011,ah012,omga001,omga002,omga011,omga012,Pa01,Pa02,ah1,ah2,a2,v,omga1,omga2,va1,vb1};
*m0,ICER,V_ICER,va,vb,Costs,dTime,cc,Pa01,Pa02,avg,a1,a2,v,ah001,ah002,va1,vb1};

create yll.&outdata from ETS[colname=ccc];
append from ETS;

quit;
%mend;
















