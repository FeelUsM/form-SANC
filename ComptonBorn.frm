*                                   p1+p2 = p3+p4
*             s - channel                                    t - channel
*
*       photon          particle                     photon    particle      
*        iu,p1           fu,p4                        iu,p1     fu,p4        
*          mu             ii                           nu        jj      
*           \             /                              \       /           
*            -           /                                -     /            
*             \   Q     /                                  \   /             
*     vert(mu,ii)-----vert(nu,ii)                      vert(nu,jj)
*             /         \                                    |
*            /           -                                   | P
*           /             \                                  |
*         ii              nu                           vert(mu,jj)
*        id,p2           fd,p3                             /   \             
*       particle        photon                            /     -            
*                                                        /       \           
*                                                       jj       mu        
*                                                     id,p2     fd,p3        
*                                                    particle  photon        

* id = vb = fu

* можно повторно делать те же самые определения
* но на множества - ругается
#redefine MyDeclareH "1"
#redefine FeynmanH "1"

#define PeskinNaumov "1"
* else SANC convention

#include myDeclare.h

#define iu  "1"
#define id  "12"
#define fu  "12"
#define fd  "1"

vector Q,P;

#if `PeskinNaumov'
GLOBAL	Born`iu'`id'`fu'`fd' =
   + Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + Ub(jj,p4,h4)*vert(1,`id',-`id',nu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',mu,jj)*U(jj,p2,h2)
	*pV(1,nu,p1,h1)*pVc(1,mu,p3,h3)
;
#else
GLOBAL	Born`iu'`id'`fu'`fd' =
   + i_*Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + i_*Ub(jj,p4,h4)*vert(1,`id',-`id',nu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',mu,jj)*U(jj,p2,h2)
	*pV(1,nu,p1,h1)*pVc(1,mu,p3,h3)
;
#endif

#call Feynman(0); * QED

* в знаменателях P и Q оставляем как есть, а вот гаммы-дирака раскрываем
id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,P) = gd(jj,p2)-gd(jj,p3);

print +s;
.sort :Compton-after-Feynman;

* для Convert()
 #define Iin1 "`iu'";   
 #define Iin2 "`id'";   
 #define Ifn3 "`fd'";   
 #define Ifn4 "`fu'";     

#$error = 0;                 * Error indicator ...*
#$num = 0;                   * Number of used momenta  (variable!) ...*
#call DiracEquation(Born`iu'`id'`fu'`fd',1) * попутно конвертируем образовавшиеся массы

.sort :Compton-after-Dirac;

sym pi,alpha,Qs,Ts;

id e^2=4*pi*alpha;
id den(1,0,Q)=1/Qs;
id den(1,0,P)=1/Ts;

print +s;
.sort :Compton-alpha;
sym volum;

#call MakeAmpSquare(Born`iu'`id'`fu'`fd',amplitudeSquared,1/4*volum)
print +s;
.sort :Compton-squred;

#$isHelicityUsed = 0;
#call Trace(amplitudeSquared,[born`iu'`id'`fu'`fd'],1,1)
print +s;
.sort :Compton-traced;

drop amplitudeSquared;
.sort :drop amplitudeSquared;;
sym Us,s,cos;

id volum=1/32/pi/s*d(cos);    
id p1.p2=1/2*(Qs+pm(`iu')^2+pm(`id')^2);
id p3.p4=1/2*(Qs+pm(`fd')^2+pm(`fu')^2);
id p2.p3=-1/2*(Ts+pm(`id')^2+pm(`fd')^2);
id p1.p4=-1/2*(Ts+pm(`iu')^2+pm(`fu')^2);
id p2.p4=-1/2*(Us+pm(`id')^2+pm(`fu')^2);
id p1.p3=-1/2*(Us+pm(`iu')^2+pm(`fd')^2);
id qel = -1;
id mgm = 0;
bracket d,den,pi,alpha,s;
print +s;
.sort :kinematics applied;

#if `PeskinNaumov'
#write "=== WE USE PESKIN-NAUMOV NOTATION ==="
#else
#write "=== WE USE SANC NOTATION ==="
#endif
.end :Compton-end;




