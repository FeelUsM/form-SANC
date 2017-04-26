*                                   p1+p2 = p3+p4
*             s - channel                                    t - channel
*
*       photon          particle                     photon    particle      
*        iu,p1           fu,p4                        iu,p1     fu,p4        
*        ii,nu           jj,mu                        ii,nu     jj,mu        
*           \             /                              \       /           
*            -           /                                -     /            
*             \   Q     /                                  \   /             
*  vert(vb,id,iu)------(vert(vb,fu,fd)               vert(vb,fu,iu)
*             /         \                                    |
*            /           -                                   | P
*           /             \                                  |
*        ii,nu           jj,mu                       vert(vb,id,fd)
*        id,p2           fd,p3                             /   \             
*       particle        photon                            /     -            
*                                                        /       \           
*                                                     ii,nu     jj,mu        
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
	*pV(1,nu,p1,h1)*pVc(1,mu,p1,h1)
;
#else
GLOBAL	Born`iu'`id'`fu'`fd' =
   + i_*Ub(ii,p4,h4)*vert(1,`id',-`id',nu,ii)*pr(`fu',Q,ii)*vert(1,`fu',-`fu',mu,ii)*U(ii,p2,h2)
	*pV(1,mu,p1,h1)*pVc(1,nu,p3,h3)
   + i_*Ub(jj,p4,h4)*vert(1,`id',-`id',nu,jj)*pr(`fu',P,jj)*vert(1,`fu',-`fu',mu,jj)*U(jj,p2,h2)
	*pV(1,nu,p1,h1)*pVc(1,mu,p1,h1)
;
#endif

#call Feynman(0); * QED

* в знаменателях P и Q оставляем как есть, а вот гаммы-дирака раскрываем
id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,P) = gd(jj,p2)-gd(jj,p3);

*#call mya2b(gd7,gd6)
*#call mya2b(g,e)

print +s;
.sort :Compton-after-Feynman;

* для Convert()
 #define Iin1 "`iu'";   * anti_fermion |
 #define Iin2 "`id'";   *      fermion | incoming
 #define Ifn3 "`fd'";   *      fermion |
 #define Ifn4 "`fu'";   * anti_fermion | outgoing  

#$error = 0;                 * Error indicator ...*
#$num = 0;                   * Number of used momenta  (variable!) ...*
#call DiracEquation(Born`iu'`id'`fu'`fd',1) 
*#call Convert(1)

.sort :Compton-after-Dirac;

sym pi,alpha,Qs,Ts;

id e^2*den(1,0,Q)=4*pi*alpha/Qs;
id e^2*den(1,0,P)=4*pi*alpha/Ts;
*id g^2/ctw^2*den(1,mp?,Q)=4*pi*alpha/Qs*4*chi(mp^2,s);
*id g^2/ctw^2*den(1,mp?,P)=4*pi*alpha/Ts*4*chi(mp^2,t);

print +s;
.sort :Compton-alpha;

#call MakeAmpSquare(Born`iu'`id'`fu'`fd',amplitudeSquared,1/4*volum)

.end :Compton-end;




