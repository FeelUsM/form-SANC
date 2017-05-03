#procedure Globals()
#ifndef `myGlobalH'
#define myGlobalH "1"

.sort :Globals-start;
#include Parameters.h

* Filling of arrays				* Globals.prc
* === Electric charges ===
Table c(1:24);
Fill c(1) = 0;
Fill c(2) = 0;
Fill c(3) = 1;
Fill c(4) = 0;
Fill c(5) = 0;
Fill c(6) = 1;
Fill c(7) = 1;
Fill c(8) =-1;
Fill c(9) = 0;
Fill c(10)= 0;
Fill c(11)= 0;
Fill c(12)= qel;
Fill c(13)= qup;
Fill c(14)= qdn;
Fill c(15)= 0;
Fill c(16)= qmo;
Fill c(17)= qch;
Fill c(18)= qst;
Fill c(19)= 0;
Fill c(20)= qta;
Fill c(21)= qtp;
Fill c(22)= qbt;
Fill c(23)= 0;
Fill c(24)= 0;

* === table of (p)article (n)ames ===
sym	w,z,h,en,el,up,dn,mn,mo,ch,st,tn,ta,tp,bt,gm,gn;			* Declar.h:88
 Table relax pn(-24:24);
  Fill pn(1) = gm;   * -- gamma
  Fill pn(-1) = gm;
  Fill pn(2) = z;    * -- Z boson
  Fill pn(-2) = z;
  Fill pn(3) = w;    * -- W bozon
  Fill pn(-3) = w;
  Fill pn(4) = h;    * -- Higgs
  Fill pn(-4) = h;
  Fill pn(11) = en;   * -- neutrino_electron
  Fill pn(-11) = en;
  Fill pn(12) = el;  * -- electron
  Fill pn(-12) = el;
  Fill pn(13) = up;  * -- up quark
  Fill pn(-13) = up;
  Fill pn(14) = dn;  * -- down quark
  Fill pn(-14) = dn;
  Fill pn(15) = mn;  * -- neutrino_muon
  Fill pn(-15) = mn;
  Fill pn(16) = mo;  * -- muon 
  Fill pn(-16) = mo;
  Fill pn(17) = ch;  * -- charmion
  Fill pn(-17) = ch;
  Fill pn(18) = st;  * -- strange quark
  Fill pn(-18) = st;
  Fill pn(19) = tn;  * -- neutrino_tau
  Fill pn(-19) = tn;
  Fill pn(20) = ta;  * -- tau
  Fill pn(-20) = ta;
  Fill pn(21) = tp;  * -- top quark
  Fill pn(-21) = tp;
  Fill pn(22) = bt;  * -- bottom quark
  Fill pn(-22) = bt;
  Fill pn(23) = gn;  * -- gluon
  Fill pn(-23) = gn;

* === masses, ratios(???), betas(???); ===
sym	mgm,mgn,mw,mz,mv,mh,men,mel,mup,mdn,mmn,mmo,mch,mst,mtn,mta,mtp,mbt,mp,mw2,mw2c; 	* Declar.h:94
* === table of (p)article (m)asses ===
 Table relax pm(-24:24);
  #do i = 1,23
     #if ((`i'<=4) || (`i'>=11))
         #$tmp = pn(`i');
         Fill pm(`i') = m`$tmp';
         Fill pm(-`i') = m`$tmp';
     #endif
  #enddo

* === кто кварк а кто лептон ===
sym	cfl,cfq;			* Declar.h:81
 Table relax clf(11:22);
  Fill clf(11) = cfl;
  Fill clf(12) = cfl;
  Fill clf(13) = cfq;
  Fill clf(14) = cfq;
  Fill clf(15) = cfl;
  Fill clf(16) = cfl;
  Fill clf(17) = cfq;
  Fill clf(18) = cfq;
  Fill clf(19) = cfl;
  Fill clf(20) = cfl;
  Fill clf(21) = cfq;
  Fill clf(22) = cfq;

 
* === tables of coupling constants ===*
* === коэффициенты правости ===
 Table vpa(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill vpa(`i') = vpa`$tmp';
  #enddo

* === коэффициенты левости ===
 Table vma(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill vma(`i') = vma`$tmp';
  #enddo



#endif
#endprocedure
