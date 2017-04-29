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

#define PeskinNaumov "1"
* else SANC convention
#$mom2mas = 1;

#if `PeskinNaumov'
#write "=== WE USE PESKIN-NAUMOV NOTATION ==="
#else
#write "=== WE USE SANC NOTATION ==="
#endif

#include Declare.h

#define iu  "1"
#define id  "12"
#define fu  "12"
#define fd  "1"

vector Q,P;

#if `PeskinNaumov'
GLOBAL	Born`iu'`id'`fu'`fd' =
   + Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + Ub(jj,p4,h4)*vert(1,`id',-`id',mu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',nu,jj)*U(jj,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
;
#else
GLOBAL	Born`iu'`id'`fu'`fd' =
   + i_*Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + i_*Ub(jj,p4,h4)*vert(1,`id',-`id',mu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',nu,jj)*U(jj,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
;
#endif

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
#call DiracEquation(Born`iu'`id'`fu'`fd',1) * попутно конвертируем образовавшиеся массы

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

bracket volum,e,mel;
print +s;
.sort :drop amplitudeSquared;;
#write "=== TRACED ==="

sym pi,alpha,Qs,Ts,Us,cos,s;

*id e^2=4*pi*alpha;
*id volum=1/32/pi/s*d(cos);    

*id den(1,m?,Q)=1/(Qs-m^2);
*id den(1,m?,P)=1/(Ts-m^2);

id p1.p2=1/2*(Qs-pm(`iu')^2-pm(`id')^2);
id p3.p4=1/2*(Qs-pm(`fd')^2-pm(`fu')^2);
id p2.p3=-1/2*(Ts-pm(`id')^2-pm(`fd')^2);
id p1.p4=-1/2*(Ts-pm(`iu')^2-pm(`fu')^2);
id p2.p4=-1/2*(Us-pm(`id')^2-pm(`fu')^2);
id p1.p3=-1/2*(Us-pm(`iu')^2-pm(`fd')^2);
id Ts = pm(`iu')^2+pm(`id')^2+pm(`fu')^2+pm(`fd')^2-Qs-Us;
id mgm = 0;

id den(1,mel,Q)=1/2/p1.p2;
id den(1,mel,P)=1/2/p3.p2;
id Qs = mel^2+2*p1.p2;
id Us = mel^2-2*p2.p3;


bracket volum,e,mel;
print +s;
.sort :kinematics applied;
#write "=== KINEMATICS ==="

#if `PeskinNaumov'
#write "=== WE USE PESKIN-NAUMOV NOTATION ==="
#else
#write "=== WE USE SANC NOTATION ==="
#endif
.end :Compton-end;




