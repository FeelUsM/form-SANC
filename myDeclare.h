#if `MyDeclareH'
#redefine MyDeclareH "0"

#define nMomenta "9";
#define nIndex "15";

sym n;
dimension n; * после этого не можем пользоваться встроенными g5_, g6_, g7_

* === импульсы ===
vector    p1,...,p`nMomenta',p;
Set sMom: p1,...,p`nMomenta',p;

* === массы ===
sym	      m1,...,m`nMomenta';
Set sMas: m1,...,m`nMomenta';

* === спиральности ===
Sym	      h1,...,h`nMomenta';
Set	sHel: h1,...,h`nMomenta';

* === set of Lorentz indices ===
Index       al,be,ga,mu,mu1,...,mu`nIndex',nu,nu1,...,nu`nIndex';
Set	sIndex: al,be,ga,mu,mu1,...,mu`nIndex',nu,nu1,...,nu`nIndex';

* === set of Spin Line Indices ===
Index     ii,jj,ii1,...,ii`nIndex';
Set sSLI: ii,jj,ii1,...,ii`nIndex';

* === gluon indices ===
dimension 8;
* attention! cgi1... reserved indices
index     gi1,...,gi`nIndex',cgi1,...,cgi`nIndex';
set sgi:  gi1,...,gi`nIndex';
set scgi:                    cgi1,...,cgi`nIndex';
set sAgi: gi1,...,gi`nIndex',cgi1,...,cgi`nIndex';
dimension n;

* === quark color indices (fundamental representation of SUc(3)) ===
dimension 3;
index     cl1,...,cl`nIndex',ccl1,...,ccl`nIndex';
set scl:  cl1,...,cl`nIndex';
set sccl:                    ccl1,...,ccl`nIndex'; * attention! ccl1... reserved indices (используется 1 раз в makeAmpSquare)
set sAcl: cl1,...,cl`nIndex',ccl1,...,ccl`nIndex';
dimension n;

*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             1 = gamma            2 = z           +/-3 = w^{+/-}    *
*====================================================================*
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

* а где множество лептонов?

Set sBI:1,2,-3,3,4,23;                                             * (s)et of (B)oson (I)ndices...
#$sBILength = 6;

Set sSI:4;                                                         * (s)et of (S)calar particle (I)ndices..
#$sSILength = 1;                                  

Set sMLBI:1,23;                                                    * (s)et of (M)ass(L)ess boson (I)ndices..
#$sMLBILength = 2;

Set sMIBI:2,-3,3;                                                  * (s)et of (M)ass(I)ve boson (I)ndices..
#$sMIBILength = 3;

* === вектора поляризации бозонов ===
* usge: pV(#bozon, mu    , momentum, helicity)
* usge: pV(sBI   , sIndex, sMom    , sHel    )
* usge: pV(gluon , mu    , momentum, helicity, gluon-index)
* usge: pV(23    , sIndex, sMom    , sHel    , sgi        )
CFun     pV,pVc;
Set spv: pV,pVc;

* === spinors ===
* usage: spinor(spin-index, momentum, helicity)
*      = spinor(sSLI      , sMom    , sHel    )
Nfun	U,Ub,V,Vb;
Set 	suf:Vb,Ub;
Set 	spf:V,U;

* === vertex ===
* usage:
*	vert(#bozon1,#bozon2,#bozon3 ,mu1   ,mu2   ,mu3   , p1,p2,p3)
*	vert(sBI    ,sBI    ,sBI     ,sIndex,sIndex,sIndex, ??,??,??)
*	vert(#bozon,#quark1 ,#quark2 ,mu    ,color1,color2,spin-index)
* = vert(sBI   ,sQI     ,sQI     ,sIndex,scl   ,scl   ,sSLI      )
*	vert(#bozon,#lepton1,#lepton2,mu                  ,spin-index)
* = vert(sBI   ,#lepton1,#lepton2,sIndex              ,sSLI      )
Nfun	vert;

* === propagator ===
* usage: 
*	pr(#bozon,mu    ,nu    ,momentum)
*	pr(sBI   ,sIndex,sIndex,sMom    )
*	pr(#quark              ,momentum,color1,color2,spin-index)
*	pr(sQI                 ,sMom    ,scl   ,scl   ,sSLI      )
*	pr(#lepton             ,momentum              ,spin-index)
*	pr(???                 ,sMom                  ,sSLI      )
Nfun	pr;

#endif
