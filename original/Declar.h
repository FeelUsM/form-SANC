#-
*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             1 = gamma            2 = z           +/-3 = w^{+/-}    *
*             4 = h                5 = phi0        +/-6 = phi^{+/-}  *
*             7 = Xp               8 = Xm             9 = Yz         *
*             10= Ya               23= g              24= Xg         *
*             11= nu_e      en     12= electron  el                  *
*             13= up        up     14= down      dn                  *
*             15= nu_mu     mn     16= muon      mu                  *
*             17= charm     ch     18= strange   st                  *
*             19= nu_tau    tn     20= tauon     ta                  *
*             21= top       tp     22= bottom    bt                  *
*====================================================================*
*
* Declarations of objects:
*
*-
#ifndef `redefineNMomenta' 
   #define nMomenta "9";
#else
   #define nMomenta "`redefineNMomenta'";
#endif
*-
#define nIndex "15";
s n,n1,n2,n3,n4,lp,ls,lt,I,FI;
*-
*- ####### QCD ########### QCD ######## QCD #######
*- 
* delta-function for QCD;
cf d;
* gluon color indices(adjoint representation of SUc(3))
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
dimension n;
*-
* attention! al1,al2,de1,de2 reserved indices
i al,be,la,ga,de,si,ka,ro;
i i,ii,j,jj,l,iu,id,fu,fd,fup,fdp,fo,fe;
i in1,...,in`nIndex',la1,...,la`nIndex',de1,...,de`nIndex';
i al1,...,al`nIndex',mu,mu1,...,mu`nIndex';
i nu,nu1,...,nu`nIndex';
*-
s AAA,BBB,CCC,DDD,higgs,higgsf,higgsb,higgsv;
s pi,alpha,alphas,lambda,omega,tHmu,tHla,volum,coef;
s e,g,gs,sr2,stw,ctw,GFermi,CFScheme,qW;
s kp,xvar,yvar,zvar,N;
* first  generation new
s qel,qup,qdn;
s ven,vel,vup,vdn,aen,ael,aup,adn,i3en,i3el,i3up,i3dn;
s vpaen,vmaen,vpael,vmael,vpaup,vmaup,vpadn,vmadn;
s dfen,dfel,dfup,dfdn;
* second generation new
s qmo,qch,qst; 
s vmn,vmo,vch,vst,amn,amo,ach,ast,i3mn,i3mo,i3ch,i3st;
s vpamn,vmamn,vpamo,vmamo,vpach,vmach,vpast,vmast;
s dfmn,dfmo,dfch,dfst;
* third  generation new
s qta,qtp,qbt;
s vtn,vta,vtp,vbt,atn,ata,atp,abt,i3tn,i3ta,i3tp,i3bt;
s vpatn,vmatn,vpata,vmata,vpatp,vmatp,vpabt,vmabt;
s dftn,dfta,dftp,dfbt;
*
s cfl,cfq;
s p1s,p2s,p3s,p4s,Qs,Ts,Us,Ps,ps,qs,pGs,pZs,eps;
s K1s,K2s,K3s,K4s,Qsmem,Tsmem,Usmem;
s qk22,qk4d2,qk4c2,int,dif,Tr,REAL;
s xi,xiz,xia,xig,xi2,xiz2,xia2,xi2m1,xiz2m1,xia2m1,xig2m1;
s s0,...,s`nMomenta',m,m0,...,m`nMomenta',x;
* masses, ratios, betas;
s w,z,h,en,el,up,dn,mn,mo,ch,st,tn,ta,tp,bt,gm,gn;
s betaw,betaz,betah,betaen,betael,betaup,betadn,betamn,betamo;
s betach,betast,betatn,betata,betatp,betabt,betagm,betagn;
s Ew,Ez,Eh,Een,Eel,Eup,Edn,Emn,Emo,Ech,Est,Etn,Eta,Etp,Ebt,Egm,Egn;
s thw,thz,thh,then,thel,thup,thdn,thmn,thmo,thch,thst,thtn,thta,thtp,thbt,thgm,thgn;
s phw,phz,phh,phen,phel,phup,phdn,phmn,phmo,phch,phst,phtn,phta,phtp,phbt,phgm,phgn;
s mgm,mgn,mw,mz,mv,mh,men,mel,mup,mdn,mmn,mmo,mch,mst,mtn,mta,mtp,mbt,mp,mw2,mw2c;
s renz,relz,rupz,rdnz,renw,relw,rupw,rdnw,renh,relh,ruph,rdnh;
s rmnz,rmoz,rchz,rstz,rmnw,rmow,rchw,rstw,rmnh,rmoh,rchh,rsth;
s rtnz,rtaz,rtpz,rbtz,rtnw,rtaw,rtpw,rbtw,rtnh,rtah,rtph,rbth;
s rwz,rhw,rhz;
s ww,wz,wtp;
*
s Za1Bos,ZmzZzBos,ZmwZwBos,ZmhZhBos,DeltarhoBos;
s Za1Lep,Za1Fer,ZmzZzFer,ZmwZwFer,ZmhZhFer,DeltarhoFer;
s Zz1Fer,Zw1Fer,Zh1Fer,ZmsqrFer,TadFer,TPallFer,calDzFer,calDwFer;
s PiZgFer,PiggFer,Za1Lep,PiggLep;
s [1/eps],[1/epsh],[p1s+mw^2],
s fix1,fix2,fix3,fix4;
s smp11,smp12,smp13,smp14,fix11,fix12,fix13,fix14;
s smp15,smp16,smp17,smp18,fix15,fix16,fix17,fix18;
s smp19,smp20,smp21,smp22,fix19,fix20,fix21,fix22;
s s,t,u,tmi,tpl,tmip,tplp,umi,upl,sqrs,spr,tau,[s-spr],cmi,cpl;
s pow,markzero;
*-
v q,q1,...,q`nMomenta',p,p1,...,p`nMomenta';
v P,Q,Qpr,k,D,S,qk3,qk2,qk4d,qk4c;
v K1,K2,K3,K4,pG,pZ,VpA,VpAc;
*-
cf Poly,Gamma,Gammai,Ftrick,pxyz;
cf comb3,combg3,comb4,gf,dgf,dgb3,dgb4,dgfbf,dgbfb,cft;
cf signf,sign,arg;
cf Tmi,Umi,Itmi,Iumi,VpV,VpVc;
cf IZ1,IZ2,IZ3,IZ4,IZ5,Z1,Z2,Z3,Z4,Z5;
cf sqrl,beta,eta,Isqrl,Ibeta,Ieta;
set sIZ:IZ1,IZ2,IZ3,IZ4,IZ5;
*-                                  REDUCTION
cf den,denc,den0,den1,den2,den3d,den3c;
cf den0i,den1i,den2i,den3di,den3ci;
cf den01,den11,den21,den3d1,den3c1;        * shifts q -> q-p1
cf den012,den112,den212,den3d12,den3c12;   * shifts q -> q-p1-p2
cf den03,den13,den23;                      * shifts q -> q-p3
cf den04,den14,den24;                      * shifts q -> q-p4
cf sym12,sym13,sym14,sym23,sym24,sym34,sym43;
cf sym112,sym122,sym211,sym221;
cf sym133,sym144,sym233,sym244,sym344,sym311,sym322,sym411,sym422,sym433;
cf sym223,sym332,sym224,sym442;
cf sym123,sym124,p1d,p2d,p3d,p4d;
cf sym11,sym22,sym33,sym44;
cf sym1112,sym1113,sym1114;
cf sym2221,sym2223,sym2224;
cf sym3331,sym3332,sym3334;
cf sym4441,sym4442,sym4443;
cf sym2111,sym3111,sym4111;
cf sym1222,sym3222,sym4222;
cf sym1333,sym2333,sym4333;
cf sym1122,sym1133,sym1144;
cf sym2211,sym2233,sym2244; 
cf sym3311,sym3322,sym3344;
cf sym4411,sym4422,sym4433; 
cf sym1123,sym1124;
cf sym2213,sym2214;
cf sym3312,sym4412;
*-
cf sym11112,sym22221,sym11122,sym22211;
cf sym111d,sym222d,sym112d,sym221d,sym1dd,sym2dd;
*-
cf sym22223,sym33332,sym22233,sym33322;
cf sym333d,sym223d,sym332d,sym3dd;
*-
cf sym22224,sym44442,sym22244,sym44422;
cf sym444d,sym224d,sym442d,sym4dd;
*-
cf sym11d,sym12d,sym13d,sym14d,sym21d,sym22d,sym23d,sym24d;
cf sym31d,sym32d,sym33d,sym34d,sym41d,sym42d,sym43d,sym44d;
*-
cf symdd;   
*-
cf invp,prop,propc,propi,propic,chi,chic;  * (p^2+s1?^2*m1?^2)
*-
cf x3i,det3,det3i,sdet3,sdet3i,f1c,f2c;    * reduction of C-functions
cf x4i,det4,det4i,sdet4,sdet4i,f1d,f2d,f3d;* reduction of D-functions
cf divisor;
*-
cf R112,R213,R232,R313,R342,R334,R356; * vector, 2-nd and 3-rd rank tensors for C-functions
*-                                     * vector, 2-nd and 3-rd rank tensors for D-functions
cf R1123,R2145,R2426,R2563,R3145,R3628,R3793,R34610,R35107,R31089,R3111213;
*-                                     * 4-th and 5-th rank tensors for C-functions
cf R513,R542,R535,R564,R556,R579,R5910,R5108,R51112;
cf R413,R442,R435,R454,R468,R487;
*-                                     * 4-th rank tensors for D-functions
cf R4145,R4627,R4893,R441013,R451311,R410614,R411158,R414712,R415129,R4131415,
cf R4161920,R4191721,R4202118;
cf Rp02112,Rp02213,Rp02232;            * vector and 2-nd rank tensors for cp02-functions
*-                                     * 3-rd, 2-nd rank tensors and vector for cp0,cp1,cp2-functions
cf Rp0313,Rp0342,Rp0334,Rp0356,Rp1313,Rp1342,Rp1334,Rp1356,Rp2313,Rp2342,Rp2334,Rp2356;
cf Rp0112,Rp0213,Rp0232,Rp1213,Rp1232,Rp1112,Rp2213,Rp2232,Rp2112;
*-
cf Rp021123,Rp022145,Rp022426,Rp022563;*vector and 2-nd rank tensors for dp02-functions
*-
cf Rp03145, Rp03628, Rp03793, Rp034610,Rp035107,Rp031089,Rp03111213; *for dp0 -functions T3,T2,V
cf Rp02145, Rp02426, Rp02563, Rp01123;
cf Rp13145, Rp13628, Rp13793, Rp134610,Rp135107,Rp131089,Rp13111213; *for dp1 -functions T3,T2,V
cf Rp12145, Rp12426, Rp12563, Rp11123;
cf Rp23145, Rp23628, Rp23793, Rp234610,Rp235107,Rp231089,Rp23111213; *for dp2 -functions T3,T2,V
cf Rp22145, Rp22426, Rp22563, Rp21123;
cf Rp33145, Rp33628, Rp33793, Rp334610,Rp335107,Rp331089,Rp33111213; *for dp3 -functions T3,T2,V
cf Rp32145, Rp32426, Rp32563, Rp31123;
*-
cf dp0221,dp0222,dp0223,dp0224,dp0225,dp0226,dp0227;                *dp02{2=T,1=V,0=S}
cf dp0224p,dp0225p,dp0226p;
cf dp0211,dp0212,dp0213,dp020;
*-                                                                  *dp0{3=T,2=T,1=V}
cf dp031,dp032,dp033,dp034,dp035,dp036,dp037,dp038,dp039,dp0310,dp0311,dp0312,dp0313;
cf dp034p,dp035p,dp036p,dp037p,dp038p,dp039p,dp0310p,dp0310pp;
cf dp021,dp022,dp023,dp024,dp025,dp026,dp027,dp011,dp012,dp013,dp00;
cf dp024p,dp025p,dp026p;
*-                                                                   dp1{3=T,2=T,1=V}
cf dp131,dp132,dp133,dp134,dp135,dp136,dp137,dp138,dp139,dp1310,dp1311,dp1312,dp1313;
cf dp134p,dp135p,dp136p,dp137p,dp138p,dp139p,dp1310p,dp1310pp;
cf dp121,dp122,dp123,dp124,dp125,dp126,dp127,dp111,dp112,dp113,dp10;
cf dp124p,dp125p,dp126p;
*-                                                                   dp2{3=T,2=T,1=V}
cf dp231,dp232,dp233,dp234,dp235,dp236,dp237,dp238,dp239,dp2310,dp2311,dp2312,dp2313;
cf dp234p,dp235p,dp236p,dp237p,dp238p,dp239p,dp2310p,dp2310pp;
cf dp221,dp222,dp223,dp224,dp225,dp226,dp227,dp211,dp212,dp213,dp20;
cf dp224p,dp225p,dp226p;
*-                                                                   dp3{3=T,2=T,1=V}
cf dp331,dp332,dp333,dp334,dp335,dp336,dp337,dp338,dp339,dp3310,dp3311,dp3312,dp3313;
cf dp334p,dp335p,dp336p,dp337p,dp338p,dp339p,dp3310p,dp3310pp;
cf dp321,dp322,dp323,dp324,dp325,dp326,dp327,dp311,dp312,dp313,dp30;
cf dp324p,dp325p,dp326p;
*-                                                                   d{4=T}
cf d41,d42,d43,d44,d45,d46,d47,d48,d49,d410,d411,d412,d413,d414,d415,d416,d417,d418,d419,d420,d421;
cf d422,d44p,d45p,d46p,d47p,d48p,d49p,d410p,d411p,d412p,d413p,d413pp,d414p,d414pp,d415p,d415pp;
cf d419p,d420p,d421p;
*-                                                                  *d{3=T}
cf d31,d32,d33,d34,d35,d36,d37,d38,d39,d310,d311,d312,d313;         
cf            d34p,d35p,d36p,d37p,d38p,d39p,d310p,d310pp;
*-                                                                  *d{2=T,1=V,0=S}
cf d21,d22,d23,d24,d25,d26,d24p,d25p,d26p,d27,d11,d12,d13,d0,ja,jaa,jav,javsub,javsubsub; 
*-
cf cp0231,cp0232,cp0233,cp0234,cp0235,cp0236; * primes              *cp02{3=T}
cf cp0221,cp0222,cp0223,cp0223p,cp0224,cp0211,cp0212,cp020;         *cp02{2=T,1=V,0=S}
cf cp031,cp032,cp033,cp034,cp035,cp036,cp033p,cp034p;               *cp0{3=T,2=T,1=V,0=S}
cf cp021,cp022,cp023,cp023p,cp024,cp011,cp012,cp00;                 
cf cp131,cp132,cp133,cp134,cp135,cp136,cp133p,cp134p;               *cp1{2=T,1=V,0=S}
cf cp121,cp122,cp123,cp123p,cp124,cp111,cp112,cp10;                 
cf cp231,cp232,cp233,cp234,cp235,cp236,cp233p,cp234p;               *cp2{2=T,1=V,0=S}
cf cp221,cp222,cp223,cp223p,cp224,cp211,cp212,cp20;                 
*-
cf c51,c52,c53,c54,c55,c56,c57,c58,c59,c510,c511,c512;              *c{5=T}
cf c53p,c54p,c55p,c56p,c59p,c510p;
cf c41,c42,c43,c44,c45,c46,c47,c48,c49,c43p,c44p,c45p,c48p;         *c{4=T}
cf c31,c32,c33,c34,c35,c36,c33p,c34p;                               *c{3=T}
cf c21,c22,c23,c24,c23p,c11,c12,c0;                                 *c{2=T,1=V,0=S}
*-
cf b43,b42,b41,b31,b32,b21,b22,b1,bff,bf,b0,b0p;
cf bpp32,bpp31,bpp22,bpp21,bpp1,bpp0;
cf bp31,bp32,bp21,bp22,bp1,bp0,bp1p,bp0p;
cf ap0,a2,a0;
*-
cf a0f,b0f,c0f,c0w,d0f,val;
*-
cf OperV,OperS,OperH;
cf Xg,Xg5,Xg6,Xg7,Yg,Yg5,Yg6,Yg7;
*-
cf qf;
*-
cf Lnmu,Lnla,ln,Li2,Re,Im,Db0f,db0f,DWb0f,dWb0f,b0d1,b0d2,b0dw,Jb,JbW,ro1ro;
cf [Qs*det3i],Invf;
*-
f gd,gd5,gd6,gd7;
f tlo,tro,tle,tre,tloc,troc,tlec,trec,U,Uc;
f vert,proj;

* -- ####################################################
* -- # Declaration of more Symbols, Functions, Sets ... #
* -- ####################################################
AutoDeclare Cf FF;                                                 * Form Factors... 
I  smIga,smIal,smIbe,smI,smImu,smInu,smImu1,...,smImu`nIndex',smInu1,...,smInu`nIndex';
I  ii1,...,ii`nIndex'; 

S  sss,sss1,...,sss9,mmm,mmm1,...,mmm9,uuu,uuu1,...,uuu9;
S  ppp,ppp1,...,ppp9,ggg,ggg1,...,ggg9,hhh,hhh1,...,hhh9;          * Symbols for internal use...
S  h1,...,h`nMomenta',four;                                        * helicities...
S  beta1,...,beta`nMomenta',E,E1,...,E`nMomenta';                  * velocity and energy...
S  ph,ph1,...,ph`nMomenta',th,th1,...,th`nMomenta';                * angles...
S  Cf,Ca,Nc;

CF pV,pVc,iZ,Z,sin,cos,sqr,iSqr,sqrt,iSqrt,sqrtL,iSqrtL,pxx,iPxx;  * various Cfunctions...
CF Nvu,Nuu,Nuv,Nvv,Npvu,Npuu,Npuv,Npvv,Nmvu,Nmuu,Nmuv,Nmvv,mV;
CF Hamp,iSqrtPx,sqrtPx,px,axpVc,axpV,sf,asf(a),dfun;               * dfun,sf,asf for internal use..
CF spinIndexList,Lambda;                                           * spinIndexList - An auxilary function for storing of spin indices
CF Tf,Ta,ff,T;
F  U,Ub,V,Vb,axgd,gd4,gdI,cR,pr;                                   * gd4 - Dirac gamma_{4}...

V  a11,...,a1`nMomenta',a21,...,a2`nMomenta',a31,...,a3`nMomenta'; * spin-vectors...
V  k0,k1;                                                          * Kleiss&Stirling vectors... 

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

Set sGP:xia,xiz,xi,xig;                                            * (s)et of (G)auge (P)aremeters...

Set sgd:gd,gd5,gd6,gd7;                   * set of Dirac Gamma matrices... 
Set sbg:g_,g5_,g6_,g7_;                   * set of Built-in Dirac Gamma matrices... 

Set suf:Vb,Ub;                            * set of anti-spinors... 
Set spf:V,U;                              * set of spinors...    

Set spv:pV,pVc;                           * set of polar. vectors...                                                
              
#ifndef `redefineIMomentaSet'
   Set sIMom:q1,...,q`nMomenta';          * set of propagator Momemta...       
#endif
#ifndef `redefineMomentaSet'
   Set sMom:p1,...,p`nMomenta',p;         * set of Momemta...       
#endif
#ifndef `redefineMassSet'
   Set sMas:m1,...,m`nMomenta';           * set of Masses...       
#endif
#ifndef `redefineHelicitySet'
   Set sHel:h1,...,h`nMomenta';           * set of helicity indices...  
#endif
#ifndef `redefineBetaSet'
   Set sBet:beta1,...,beta`nMomenta';     * set of betas... 
#endif
#ifndef `redefineEnergySet'
   Set sEng:E1,...,E`nMomenta';           * set of Energies...  
#endif
#ifndef `redefineAngleSet'
   Set sPh:ph1,...,ph`nMomenta';          * set of ph's...  
   Set sTh:th1,...,th`nMomenta';          * set of th's...
#endif

Set sA1:a11,...,a1`nMomenta';             * 
Set sA2:a21,...,a2`nMomenta';             * set of spin-vectors...
Set sA3:a31,...,a3`nMomenta';             * 
Set sASV:a11,...,a1`nMomenta',
         a21,...,a2`nMomenta',
         a31,...,a3`nMomenta';

Set sIndex:ga,al,be,mu,mu1,...,mu`nIndex',nu,nu1,...,nu`nIndex';             * set of Lorentz indices...
Set sIndexSm:smIga,smIal,smIbe,smImu,smImu1,...,smImu`nIndex',smInu,smInu1,...,smInu`nIndex';  * ...
Set sSLI:ii,jj,ii1,...,ii`nIndex';              * set of Spin Line Indices...

#$idDotProductsUsed = 0;     * Check whether "idDotProducts.prc" is used ...*
#$isHelicityUsed = 0;        * Check whether helicity is used ...*
#$num = 0;                   * Number of used momenta  (variable!) ...*
#$error = 0;                 * Error indicator ...*
#$RSMomentumPair = -1;       * Rest System Momentum Pair (see TraceHelicity.prc)...* 
#define firstTime "0"; 
*-
nwrite statistics;
.global
#+
