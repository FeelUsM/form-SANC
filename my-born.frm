* === общие символы ===
sym	n,n1,n2,n3,n4,lp,ls,lt,I,FI;			* Decalr.h:23

* === индексы ===
#define nIndex "15";
index	cl1,...,cl`nIndex',
	ccl1,...,ccl`nIndex';				* Decalr.h:45
index	al,be,la,ga,de,si,ka,ro;			* Decalr.h:55
index	i,ii,j,jj,l,iu,id,fu,fd,fup,fdp,fo,fe;		* Decalr.h:56
index	in1,...,in`nIndex',
	la1,...,la`nIndex',
	de1,...,de`nIndex';
index	al1,...,al`nIndex',
	mu,mu1,...,mu`nIndex',
	nu,nu1,...,nu`nIndex';				* Decalr.h:59
Index	ii1,...,ii`nIndex'; 				* Declar.h:269
Set 	sSLI:ii,jj,ii1,...,ii`nIndex';              * set of Spin Line Indices...
Set	sIndex:ga,al,be,mu,mu1,...,mu`nIndex',nu,nu1,...,nu`nIndex';             * set of Lorentz indices...
Index	smIga,smIal,smIbe,smI,smImu,smInu,smImu1,...,smImu`nIndex',smInu1,...,smInu`nIndex';
Set	sIndexSm:smIga,smIal,smIbe,smImu,smImu1,...,smImu`nIndex',smInu,smInu1,...,smInu`nIndex';  * ...

* === векторы ===
#define nMomenta "9";
vector	q,q1,...,q`nMomenta',p,p1,...,p`nMomenta';	* Declar.h:113
Set 	sMom:p1,...,p`nMomenta',p;         * set of Momemta...       
vector	P,Q,Qpr,k,D,S,qk3,qk2,qk4d,qk4c;

* === спиральности ===
Sym	h1,...,h`nMomenta',four; 			* Declar.h:273
Set	sHel:h1,...,h`nMomenta';

CF pV,pVc,iZ,Z,sin,cos,sqr,iSqr,sqrt,iSqrt,sqrtL,iSqrtL,pxx,iPxx;  * various Cfunctions...
Set	spv:pV,pVc;	* вектора поляризации бозонов

* === Дираковские матрицы и спиноры, и то, что их может содержать ===
Nfun	axgd,gd4,gdI,cR,pr;				* Declar.h:283          * gd4 - Dirac gamma_{4}...
Nfun	gd,gd5,gd6,gd7;					* Declar.h:260
Set sgd:gd,gd5,gd6,gd7;                   * set of Dirac Gamma matrices... 
Set sbg:g_,g5_,g6_,g7_;                   * set of Built-in Dirac Gamma matrices... 

Nfun	U,Ub,V,Vb,vert,proj;				* Declar.h:262
Set 	suf:Vb,Ub;                            * set of anti-spinors... 
Set 	spf:V,U;                              * set of spinors...    

#define iu  "12"
#define id  "12"
#define fu  "12"
#define fd  "12"

* === Поехали ===

* vb={1,2} - gamma, z
* ii,jj,ii1,ii2 - спиновые индексы

GLOBAL Born`iu'`id'`fu'`fd' =
#do vb={1,2} 
   + i_*Vb(ii,p1,h1)*vert(-`vb',`id',-`iu',nu,ii)*U(ii,p2,h2)*pr(`vb',nu,mu,Q)
       *Ub(jj,p3,h3)*vert(`vb',`fu',-`fd',mu,jj)*V(jj,p4,h4)
   - i_*Ub(ii2,p3,h3)*vert(-`vb',`id',-`fd',mu1,ii2)*U(ii2,p2,h2)*pr(`vb',nu1,mu1,P)
       *Vb(ii1,p1,h1)*vert(`vb',`fu',-`iu',nu1,ii1)*V(ii1,p4,h4)
#enddo
;
#write "=== born created: ==="
print +s;
.sort:born created;
#write "=== -born created: ==="

* ###########################################################################################
* #                                                                                         #
* #                                FEYNMAN RULES FOR SM                                     #
* #                                                                                         #
* ###########################################################################################

sym	e,g,gs,sr2,stw,ctw,GFermi,CFScheme,qW;

* === delta-function for QCD - CKM-матрица ===
cfun	d;					* Declar.h:31	

dimension 8;
* attention! cgi1... reserved indices
i gi1,...,gi`nIndex',cgi1,...,cgi`nIndex';
set sgi: gi1,...,gi`nIndex';
set scgi: cgi1,...,cgi`nIndex';
set sAgi:gi1,...,gi`nIndex',cgi1,...,cgi`nIndex';
*
*-
* quark color indices(fundamental representation of SUc(3))
dimension 3;
* attention! ccl1... reserved indices
i cl1,...,cl`nIndex',ccl1,...,ccl`nIndex';
set scl: cl1,...,cl`nIndex';
set sccl:ccl1,...,ccl`nIndex';
set sAcl: cl1,...,cl`nIndex',ccl1,...,ccl`nIndex';
*-
*- ####### QCD ########### QCD ######## QCD #######
*-
dimension n; * после этого не можем пользоваться встроенными g5_, g6_, g7_
*-

* === first  generation new ===			* Declar.h:65
sym	qel,qup,qdn;
sym	ven,vel,vup,vdn,aen,ael,aup,adn,i3en,i3el,i3up,i3dn;
sym	vpaen,vmaen,vpael,vmael,vpaup,vmaup,vpadn,vmadn;
sym	dfen,dfel,dfup,dfdn;
* === second generation new ===
sym	qmo,qch,qst; 
sym	vmn,vmo,vch,vst,amn,amo,ach,ast,i3mn,i3mo,i3ch,i3st;
sym	vpamn,vmamn,vpamo,vmamo,vpach,vmach,vpast,vmast;
sym	dfmn,dfmo,dfch,dfst;
* === third  generation new ===
sym	qta,qtp,qbt;
sym	vtn,vta,vtp,vbt,atn,ata,atp,abt,i3tn,i3ta,i3tp,i3bt;
sym	vpatn,vmatn,vpata,vmata,vpatp,vmatp,vpabt,vmabt;
sym	dftn,dfta,dftp,dfbt;

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

* === masses, ratios, betas; ===
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

* === sets of indices ===
Set sFI:-11,11,-12,12,-13,13,-14,14,-15,15,-16,16,                 
        -17,17,-18,18,-19,19,-20,20,-21,21,-22,22;                 * (s)et of (F)ermion (I)ndices...
#$sFILength = 24;
Set sIFI: 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22;          * (s)et of (Incoming) (F)ermion (I)ndices...
#$sIFILength = 12;
Set sOFI:-11,-12,-13,-14,-15,-16,-17,-18,-19,-20,-21,-22;          * (s)et of (Outgoing) (F)ermion (I)ndices...
#$sOFILength = 12;

Set sQI:-13,13,-14,14,-17,17,-18,18,-21,21,-22,22;                 * (s)et of (Q)uark (I)ndices...
#$sQILength = 11;
Set sIQI: 13, 14, 17, 18, 21, 22;                                  * (s)et of (Incoming) (Q)uark (I)ndices...
#$sIQILength = 6;
Set sOQI:-13,-14,-17,-18,-21,-22;                                  * (s)et of (Outgoing) (Q)uark (I)ndices...
#$sOQILength = 6;

Set sBI:1,2,-3,3,4,23;                                             * (s)et of (B)oson (I)ndices...
#$sBILength = 6;

Set sSI:4;                                                         * (s)et of (S)calar particle (I)ndices..
#$sSILength = 1;                                  

Set sMLBI:1,23;                                                    * (s)et of (M)ass(L)ess boson (I)ndices..
#$sMLBILength = 2;

Set sMIBI:2,-3,3;                                                  * (s)et of (M)ass(I)ve boson (I)ndices..
#$sMIBILength = 3;

* === Поехали ===

*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             1 = gamma            2 = z           +/-3 = w^{+/-}    *
*====================================================================*
*             -------------boson_boson_vertices-----------
* === gamma,w+,w- ===
id vert(n?{1,-1},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{1,-1},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
* === z,w+,w- ===
id vert(n?{2,-2},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{2,-2},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
*
* ((1,2)33)--(pqk) ---> (kpq)--(33(1,2))
id vert(+3,-3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,+3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
* ((1,2)33)--(pqk) ---> (qkp)--(3(1,2)3)
id vert(+3,n?{1,-1},-3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))    
                                                          +d_(mu,be)*(k(al)-p(al))); 
id vert(-3,n?{1,-1},+3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))      
                                                          +d_(mu,be)*(k(al)-p(al))); 
id vert(+3,-3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =    +g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))    
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,+3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =    -g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))     
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,n?{2,-2},+3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
id vert(+3,n?{2,-2},-3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             11= nu_e      en     12= electron  el                  *
*             13= up        up     14= down      dn                  *
*             15= nu_mu     mn     16= muon      mu                  *
*             17= charm     ch     18= strange   st                  *
*             19= nu_tau    tn     20= tauon     ta                  *
*             21= top       tp     22= bottom    bt                  *
*====================================================================*
*
*-             ----------boson_fermion_vertices-----------
*
* first generation
* ----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+11,-11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-11,+11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+12,-12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-12,+12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+13,-13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-13,+13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+14,-14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-14,+14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{2,-2},+11,-11,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},-11,+11,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},+12,-12,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},-12,+12,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},+13,-13,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},-13,+13,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},+14,-14,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
id vert(n?{2,-2},-14,+14,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
*
id vert(-3,+11,-12,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-12,+11,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+12,-11,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-11,+12,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+13,-14,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-14,+13,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+14,-13,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-13,+14,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
* second generation
* -----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+15,-15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-15,+15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+16,-16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-16,+16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+17,-17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-17,+17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+18,-18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-18,+18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{2,-2},+15,-15,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},-15,+15,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},+16,-16,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},-16,+16,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},+17,-17,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},-17,+17,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},+18,-18,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
id vert(n?{2,-2},-18,+18,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
*
id vert(-3,+15,-16,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-16,+15,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+16,-15,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-15,+16,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+17,-18,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-18,+17,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+18,-17,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-17,+18,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
* third generation
* ----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n1)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n2)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4      /ctw*gd(ii,mu)*(vpa(n1)*gd6(ii)+vma(n1)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4      /ctw*gd(ii,mu)*(vpa(n2)*gd6(ii)+vma(n2)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+21    ,-22    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-22    ,+21    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+22    ,-21    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-21    ,+22    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
* SM vertices
id vert(n?{1,-1},+19,-19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-19,+19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+20,-20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-20,+20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+21,-21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-21,+21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+22,-22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-22,+22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{2,-2},+19,-19,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},-19,+19,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},+20,-20,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},-20,+20,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},+21,-21,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},-21,+21,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},+22,-22,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
id vert(n?{2,-2},-22,+22,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
*
id vert(-3,+19,-20,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-20,+19,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+20,-19,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-19,+20,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+21,-22,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-22,+21,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+22,-21,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-21,+22,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
.sort:vertices expanded;
*-                                  REDUCTION
cfun	den,denc,den0,den1,den2,den3d,den3c;				* Declar.h:125
* (s)et of (G)auge (P)aremeters...
sym	xi,xiz,xia,xig,xi2,xiz2,xia2,xi2m1,xiz2m1,xia2m1,xig2m1;	* Declar.h:85
Set	sGP:xia,xiz,xi,xig;           					* Declar.h:315
* table of gauge parameters *
 Table relax rxi(1:23);
  Fill rxi(1) = xia;
  Fill rxi(2) = xiz;
  Fill rxi(3) = xi;
  Fill rxi(23)= xig;

* ###################################################
* #            P R O P A G A T O R S                #
* ###################################################
*
*   
* ========== SM Bosons ==========
* ===============================
*
*id pr(n?{1,-1},mu?,nu?,p?)  = den(1,0,p)*(d_(mu,nu) + (rxi(1)^2-1)*p(mu)*p(nu)*den(1,0,p));
*id pr(n?{2,-2},mu?,nu?,p?)  = den(1,pm(2),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(2)^-1*pm(2)^-1)
*                              - p(mu)*p(nu)*pm(2)^-1*pm(2)^-1*den(rxi(2),pm(2),p);
*id pr(n?{3,-3},mu?,nu?,p?)  = den(1,pm(3),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(3)^-1*pm(3)^-1)
*                              - p(mu)*p(nu)*pm(3)^-1*pm(3)^-1*den(rxi(3),pm(3),p);

id pr(n?{1,-1},mu?,nu?,p?)  = den(1,0,p)*(d_(mu,nu));
id pr(n?{2,-2},mu?,nu?,p?)  = den(1,pm(2),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(2)^-1*pm(2)^-1);
id pr(n?{3,-3},mu?,nu?,p?)  = den(1,pm(3),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(3)^-1*pm(3)^-1);

 repeat;
   id d(?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = d(?a,?c,?d,?b);
   id d(?a,cl1?sAgi,?b)*d(?c,cl1?sAgi,?d) = d(?a,?c,?d,?b);
 endrepeat;

#define I "1";
* I=0: QED,QCD
* I=1: EW
#if {`I'} = 0
    #call a2b(gd6,gd5)
    #call a2b(gd7,gd5)
    #call a2b(g,e)
#endif

#write "=== FeynmanRules applied: ==="
print +s;
.sort:FeynmanRules applied;
#write "=== -FeynmanRules applied ==="

id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,Q) = gd(jj,p3)+gd(jj,p4);
id gd(ii1,P)=-gd(ii1,p1)+gd(ii1,p4);
id gd(ii2,P)= gd(ii2,p2)-gd(ii2,p3);
#write "=== Q and P expanded: ==="
print +s;
.sort:Q and P expanded;
#write "=== -Q and P expanded ==="

Cfun	sf;	* for internal use
Sym	sss,sss1,...,sss9,mmm,mmm1,...,mmm9,uuu,uuu1,...,uuu9;
Sym	ppp,ppp1,...,ppp9,ggg,ggg1,...,ggg9,hhh,hhh1,...,hhh9;          * Symbols for internal use...
sym	s0,...,s`nMomenta',m,m0,...,m`nMomenta',x;
Set sMas:m1,...,m`nMomenta';           * set of Masses...       

 #define Iin1 "`iu'";   * anti_fermion |
 #define Iin2 "`id'";   *      fermion | incoming
 #define Ifn3 "`fd'";   *      fermion |
 #define Ifn4 "`fu'";   * anti_fermion | outgoing  

#$error = 0;                 * Error indicator ...*
#call DiracEquation(Born`iu'`id'`fu'`fd',-1)
#$num = 0;                   * Number of used momenta  (variable!) ...*
#call Convert(1)

#write "=== DiracEquation applied and converted: ==="
print +s;
.sort:DiracEquation applied and converted;
#write "=== -DiracEquation applied and converted ==="

sym	p1s,p2s,p3s,p4s,Qs,Ts,Us,Ps,ps,qs,pGs,pZs,eps;			* Declar.h:82
sym	pi,alpha,alphas,lambda,omega,tHmu,tHla,volum,coef;		* Declar.h:62
cfun	invp,prop,propc,propi,propic,chi,chic;  			* Declar.h:166	* (p^2+s1?^2*m1?^2)
sym	s,t,u,tmi,tpl,tmip,tplp,umi,upl,sqrs,spr,tau,[s-spr],cmi,cpl;	* Declar.h:110

#call a2b(gd7,gd6)
id den(1,mh,Q)=prop(mh^2,Qs);
id den(1,mh,P)=prop(mh^2,Ts);
id g^2*stw^2*den(1,0,Q)=4*pi*alpha/Qs;
id g^2*stw^2*den(1,0,P)=4*pi*alpha/Ts;
id g^2/ctw^2*den(1,mp?,Q)=4*pi*alpha/Qs*4*chi(mp^2,s);
id g^2/ctw^2*den(1,mp?,P)=4*pi*alpha/Ts*4*chi(mp^2,t);
#write "=== Denumerators expanded: ==="
print +s;
.sort:Denumerators expanded;
#write "=== -Denumerators expanded ==="

CF Hamp,iSqrtPx,sqrtPx,px,axpVc,axpV,sf,asf(a),dfun;               * dfun,sf,asf for internal use..
CF Tf,Ta,ff,T;
S  Cf,Ca,Nc;

#call MakeAmpSquare(Born`iu'`id'`fu'`fd',amplitudeSquared,1/4*volum)
#write "=== Squared: ==="
print +s;
.sort:Squared;
#write "=== -Squared ==="
