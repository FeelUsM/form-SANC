*                                   p1+p2 = p3+p4
*             s - channel                                    t - channel
*
*       photon          particle                     photon    particle      
*        iu,p1           fu,p4                        iu,p1     fu,p4        
*          mu             ii                           mu        jj      
*           \             /                              \       /           
*            -           /                                -     /            
*             \   Q     /                                  \   /             
*     vert(mu,ii)-----vert(nu,ii)                      vert(mu,jj)
*             /         \                                    |
*            /           -                                   | P
*           /             \                                  |
*         ii              nu                           vert(nu,jj)
*        id,p2           fd,p3                             /   \             
*       particle        photon                            /     -            
*                                                        /       \           
*                                                       jj       nu        
*                                                     id,p2     fd,p3        
*                                                    particle  photon        

#include Declar.h
#call Globals();

#define iu  "1"
#define id  "12"
#define fu  "12"
#define fd  "1"

vector Q,P;

GLOBAL	Born`iu'`id'`fu'`fd' =
   + i_*Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + i_*Ub(jj,p4,h4)*vert(1,`id',-`id',mu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',nu,jj)*U(jj,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
;

#call FeynmanRules(0); * QED
id qel = -1;

* в знаменателях P и Q оставляем как есть, а вот гаммы-дирака раскрываем
id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,P) = gd(jj,p2)-gd(jj,p3);

bracket e;
print +s;
.sort :Compton-after-Feynman;
#write "=== FEYNMAN ==="

* для Convert()
 #define Iin1 "`iu'";   
 #define Iin2 "`id'";   
 #define Ifn3 "`fd'";   
 #define Ifn4 "`fu'";     

#$error = 0;                 * Error indicator ...*
#$num = 0;                   * Number of used momenta  (variable!) ...*
#ifndef `useDirac'
  #define useDirac "1"
#endif
#if `useDirac'
  #call DiracEquation(Born`iu'`id'`fu'`fd',1) * попутно конвертируем образовавшиеся массы
#endif

bracket e;
print +s;
.sort :Compton-after-Dirac;
#write "=== DIRAC ==="

sym volum;

#call MakeAmpSquare(Born`iu'`id'`fu'`fd',amplitudeSquared,1/4*volum)
bracket volum,e,den;
print +s;
.sort :Compton-squred;
#write "=== SQURED ==="

#$isHelicityUsed = 0;
#call Trace(amplitudeSquared,[born`iu'`id'`fu'`fd'],1,1)
.sort :Compton-traced;
drop amplitudeSquared;

bracket volum,e,den;
print +s;
.sort :drop amplitudeSquared;;
#write "=== TRACED ==="

sym pi,alpha,SS,TT,UU;

*id e^2=4*pi*alpha;
*id volum=1/32/pi/SS*d(cos);    

*id den(1,m?,Q)=1/(Qs-m^2);
*id den(1,m?,P)=1/(Ts-m^2);

id p1.p2=1/2*(SS-pm(`iu')^2-pm(`id')^2);
id p3.p4=1/2*(SS-pm(`fd')^2-pm(`fu')^2);
id p2.p3=-1/2*(UU-pm(`id')^2-pm(`fd')^2);
id p1.p4=-1/2*(UU-pm(`iu')^2-pm(`fu')^2);
id p2.p4=-1/2*(TT-pm(`id')^2-pm(`fu')^2);
id p1.p3=-1/2*(TT-pm(`iu')^2-pm(`fd')^2);
id TT = pm(`iu')^2+pm(`id')^2+pm(`fu')^2+pm(`fd')^2-SS-UU;
id mgm = 0;

bracket volum,e,den;
print +s;
.sort :kinematics applied-1;
#write "=== KINEMATICS-1 ==="

id den(1,mel,Q)=1/2/p1.p2;
id den(1,mel,P)=-1/2/p3.p2;
id SS = mel^2+2*p1.p2;
id UU = mel^2-2*p2.p3;

bracket volum,e,mel;
print +s;
.sort :kinematics applied-2;
#write "=== KINEMATICS-2 ==="

.end :Compton-end;




