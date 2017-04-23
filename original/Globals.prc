#procedure Globals(dummy)
*------------------------
*
* Filling of arrays
*
Table c(1:24);
*
* Electric charges
*
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
*
Table relax cd(-24:24);
*
* Electric charges for \delta-functions
*
Fill cd(0) = 0;
Fill cd(1) = 0;
Fill cd(2) = 0;
Fill cd(3) = 1;
Fill cd(4) = 0;
Fill cd(5) = 0;
Fill cd(6) = 1;
Fill cd(7) = 1;
Fill cd(8) =-1;
Fill cd(9) = 0;
Fill cd(10)= 0;
Fill cd(11)= 0;
Fill cd(12)=-1;
Fill cd(13)=+2/3;
Fill cd(14)=-1/3;
Fill cd(15)= 0;
Fill cd(16)=-1;
Fill cd(17)=+2/3;
Fill cd(18)=-1/3;
Fill cd(19)= 0;
Fill cd(20)=-1;
Fill cd(21)=+2/3;
Fill cd(22)=-1/3;
Fill cd(23)= 0;
Fill cd(24)= 0;
*
Fill cd(-1) = 0;
Fill cd(-2) = 0;
Fill cd(-3) =-1;
Fill cd(-4) = 0;
Fill cd(-5) = 0;
Fill cd(-6) =-1;
Fill cd(-7) =-1;
Fill cd(-8) = 1;
Fill cd(-9) = 0;
Fill cd(-10)= 0;
Fill cd(-11)= 0;
Fill cd(-12)= 1;
Fill cd(-13)=-2/3;
Fill cd(-14)=+1/3;
Fill cd(-15)= 0;
Fill cd(-16)= 1;
Fill cd(-17)=-2/3;
Fill cd(-18)=+1/3;
Fill cd(-19)= 0;
Fill cd(-20)= 1;
Fill cd(-21)=-2/3;
Fill cd(-22)=+1/3;
Fill cd(-23)= 0;
Fill cd(-24)= 0;
*
Table gh(0:23);
*
* ghost's charges
*
Fill gh(0) = 0;
Fill gh(1) = 0;
Fill gh(2) = 0;
Fill gh(3) = 0;
Fill gh(4) = 0;
Fill gh(5) = 0;
Fill gh(6) = 0;
Fill gh(7) = 1;
Fill gh(8) =-1;
Fill gh(9) = 1;
Fill gh(10)= 1;
Fill gh(11)= 0;
Fill gh(12)= 0;
Fill gh(13)= 0;
Fill gh(14)= 0;
Fill gh(15)= 0;
Fill gh(16)= 0;
Fill gh(17)= 0;
Fill gh(18)= 0;
Fill gh(19)= 0;
Fill gh(20)= 0;
Fill gh(21)= 0;
Fill gh(22)= 0;
Fill gh(23)= 1;

* table of (p)article (n)ames *
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

* table of (p)article (m)asses *
 Table relax pm(-24:24);
  #do i = 1,23
     #if ((`i'<=4) || (`i'>=11))
         #$tmp = pn(`i');
         Fill pm(`i') = m`$tmp';
         Fill pm(-`i') = m`$tmp';
     #endif
  #enddo

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

* table of (rat)ios *
 Table relax rat(11:22,2:4);
  #do i = {2,3,4}
     #$tmp1 = pn(`i'); 
     #do j = 11,22 
        #$tmp2 = pn(`j');        
        Fill rat(`j',`i') = r`$tmp2'`$tmp1';
     #enddo
  #enddo
 
* tables of coupling constants *
 Table vpa(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill vpa(`i') = vpa`$tmp';
  #enddo

 Table vma(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill vma(`i') = vma`$tmp';
  #enddo

 Table relax vf(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill vf(`i') = v`$tmp';
  #enddo

 Table relax af(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill af(`i') = a`$tmp';
  #enddo

 Table relax df(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill df(`i') = df`$tmp';
  #enddo

* table of weak isospin *
 Table i3(11:22);
  #do i = 11,22  
     #$tmp = pn(`i');
     Fill i3(`i') = i3`$tmp';
  #enddo

* table of gauge parameters *
 Table relax rxi(1:23);
  Fill rxi(1) = xia;
  Fill rxi(2) = xiz;
  Fill rxi(3) = xi;
  Fill rxi(23)= xig;

 #ifdef `debug'
    printtable c;
    printtable cd;
    printtable gh;
    printtable pn;
    printtable pm;
    printtable clf;
    printtable rat;
    printtable vpa;
    printtable vma;
    printtable vf;
    printtable af;
    printtable i3;
    printtable rxi;
 #endif 
 
.global
#endprocedure
*------------ 
