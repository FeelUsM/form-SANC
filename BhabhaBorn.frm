*---------------------------------
*                                   p1+p2 = p3+p4
*             s - channel                                    t - channel
*
*     anti-particle  anti-particle                anti-particle          anti-particle
*        iu,p1           fu,p4                     iu,p1,ii,nu            fu,p4,ii,nu  
*        ii,nu           jj,mu                                 vert(vb,fu,iu) 
*           \             /                              <----------/<----------
*            \           /                                          \         
*             \   Q     /                                           / P
*  vert(vb,id,iu)/\/\/\/(vert(vb,fu,fd)                             \
*             /         \                                           /
*            /           \                               ---------->\---------->
*           /             \                                   vert(vb,id,fd) 
*        ii,nu           jj,mu                           
*        id,p2           fd,p3                      id,p2,jj,mu           fd,p3,jj,mu 
*       particle       particle                      particle               particle
*
*

#ifndef `PeskinNaumov'
  #define PeskinNaumov "1"
#endif
* else SANC convention
#$mom2mas = 1;

#if `PeskinNaumov'
#write "=== WE USE PESKIN-NAUMOV NOTATION ==="
#else
#write "=== WE USE SANC NOTATION ==="
#endif

#-
#include Declare.h
#+

#define iu  "12"
#define id  "12"
#define fu  "16"
#define fd  "16"

vector Q,P;

GLOBAL Born`iu'`id'`fu'`fd' =
#do vb={1,2}
*#define vb "1"
   + Vb(ii ,p1,h1)*vert(-`vb',`id',-`iu',nu ,ii )*U(ii ,p2,h2)*pr(`vb',nu ,mu ,Q)
    *Ub(jj ,p3,h3)*vert( `vb',`fu',-`fd',mu ,jj )*V(jj ,p4,h4)
   - Ub(ii2,p3,h3)*vert(-`vb',`id',-`fd',mu1,ii2)*U(ii2,p2,h2)*pr(`vb',nu1,mu1,P)
    *Vb(ii1,p1,h1)*vert( `vb',`fu',-`iu',nu1,ii1)*V(ii1,p4,h4)
#enddo
;
#write "=== born created ==="
print +s;
.sort:born created;

#call FeynmanRules(1); * EW считет быстрее
id qel = -1;
id qmo = -1;
id vert(?e) = 0;

* в знаменателях P и Q оставляем как есть, а вот гаммы-дирака раскрываем
id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,Q) = gd(jj,p3)+gd(jj,p4);
id gd(ii1,P)=-gd(ii1,p1)+gd(ii1,p4);
id gd(ii2,P)= gd(ii2,p2)-gd(ii2,p3);

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
  #define useDirac "0"
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

sym pi,alpha,SS,TT,UU,cos;

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
id mel = 0;

bracket volum,e,den;
print +s;
.sort :kinematics applied-1;
#write "=== KINEMATICS-1 ==="

*id SS = mel^2+2*p1.p2;
id UU = mmo^2-2*p2.p3;

bracket volum,e,mmo,den;
print +s;
.sort :kinematics applied-2;
#write "=== KINEMATICS-2 ==="

#if `PeskinNaumov'
#write "=== WE USE PESKIN-NAUMOV NOTATION ==="
#else
#write "=== WE USE SANC NOTATION ==="
#endif
.end :Compton-end;




