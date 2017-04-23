*============================================================================
*                                                                           
*                File: Brem_Bhabha.frm                                  
*============================================================================
* Calculation of the Photon Bremsstrahlung for s-channel and t-channel Bhabha
*============================================================================
*                                                                        
*     antifermion        antifermion           antifermion      antifermion
*    `iu',p1,ii,mu      `fu',p4,jj,nu         `iu',p1,ii,mu    `fu',p4,jj,nu
*            \             /                         \             /
*             \           /                           \           /
*              \  B,phi  /                             \  B,phi  /
* vert(vn,fo,-iu)-------vert(vn,fu,-fd)   vert(vn,id,-iu)-------vert(vn,fu,-fe)
*              /|\ Q     \                             /|  Q    |
*               |         \                           /         |
*            fo |          jj,nu                      ii,mu     | fe 
*         -p5+p2|          `fd',p3                `id',p2       | p5+p3
*               |          fermion                fermion       |
*               |                                               |
*               |                                               |
*               |      p5                                      \|/  photon
* vert(ph,id,-fo)--------                         vert(ph,fe,-fd)---------
*              /|  photon                                       \     p5
*             /                                                  \
*            /                                                    \
*           /                                                      \
*      `id',p2,ii,al                                            `fd',p3,jj,al       
*         fermion                                                  fermion        
*
*
* 
*     antifermion                                                antifermion            
*     `iu',p1,ii,al                                              `fu',p4,jj,al
*            \                                                     /
*             \                                                   /
*              \  photon                                         /      p5
* vert(ph,fo,-iu)--------                         vert(ph,fu,-fe)---------
*              /|\    p5                                        |   photon
*               |                                               |
*               |                                               |
*               |                                               |
*               |      antifermion           antifermion        |
*               |     `fu',p4,jj,nu         `iu',p1,ii,mu       | fe
*            fo |          /                         \          | -p5-p4
*         p5-p1 |         /                           \         |
*               | B,phi  /                             \ B,phi \|/
* vert(vs,id,-fo)-------vert(vs,fu,-fd)   vert(vs,id,-iu)-------vert(vs,fe,-fd)
*              /   Q     \                             /  Q      \
*             /           \                           /           \
*            /             \                         /             \
*    `id',p2,ii,mu     `fd',p3,jj,nu          `id',p2,ii,mu      `fd',p3,jj,nu
*      fermion           fermion                fermion            fermion
*
*                                 
*  Born:            p1+p2=p3+p4
*
*  Q = p1+p2 = p3+p4          P=p2-p3=p4-p1      Us=(p2-p4)^2=(p1-p3)^2 
*  Q^2=Qs= -s;  sqrs^2=s;     P^2=Ts;
*  id volum=1/32/pi/s*d(cos);   
*  Qs=-s;   beta(-Qs,mel^2,mel^2)^2 = beta(s,mel^2,mel^2)^2= 1-4*mel^2/s; 
*  Ts=-2*mel^2+s/2*(1-beta^2*cos(theta)) -> Ts=s/2*(1-cos(theta))
*
*  Bremsstrahlung:             p1+p2=p3+p4+p5
*     Z1=-2*p1.p5=-2*p1.n*p5^0,               Z2=-2*p2.p5=-2*p2.n*p5^0,     
*     Z3=-2*p3.p5=-2*p3.n*p5^0,               Z4=-2*p4.p5=-2*p4.n*p5^0.
*    
*  e^+(p1) + e^-(p2) ==> e^-(p3) + e^-(p4) + photon(p5)
*  m1=m2=m3=m4=mel;  p1^0=p2^0=k3^0=k4^0=sqr/2;   beta=beta(s,mel^2,mel^2)
*  beta^2=1-4*mel^2/s;
*  sqrlT^2=Ts^2+4*Ts*mel^2     sqrlU^2=Us^2+4*Us*mel^2
*  
* Re(Jb(-s,m,m))= 1/2*sqrl*\int\limits_{0}^{1}dy 
*                                       \frac{1}{-s*y*(1-y)+m^2*(1-y)+m^2*y} 
* Re(JbW(Ts,m1,m2))= -sqrl*\int\limits_{0}^{1}dy 
*                                       \frac{1}{Ts*y*(1-y)+m1^2*(1-y)+m2^2*y}  
*
* pV - polarization vector of the photon
*=============================================================================
cf IZ1,IZ2,IZ3,IZ4,Z1,Z2,Z3,Z4;
#include Declar.h
#call Globals()
*
#define typeIU  "12"
#define typeID  "12"
#define typeFU  "12"
#define typeFD  "12"

#ifdef `typeIU'
#ifdef `typeID'
#ifdef `typeFU'
#ifdef `typeFD'
*
 Set sLeftTreeMomenta  : p3;
 Set sRightTreeMomenta : p4,p5;
 Set sCMSMomenta       : p2,p1;
 #define Iin1 "`typeIU'";   * anti_fermion |
 #define Iin2 "`typeID'";   *      fermion | incoming
 #define Ifn3 "`typeFD'";   *      fermion |
 #define Ifn4 "`typeFU'";   * anti_fermion | outgoing  
 #define Ifn5 "1";          * photon       | outgoing  
*
* program
*--------
.sort
s BORN`typeIU'`typeID'`typeFU'`typeFD';
#call BornBhabha(`typeIU',`typeID',`typeFU',`typeFD')
.sort
on stat;
id mel=0;
b g,stw,ctw,volum,prop,qel,chi,chic,pi,alpha,vmael,vpael;
*print Born`typeIU'`typeID'`typeFU'`typeFD';
*print [born`typeIU'`typeID'`typeFU'`typeFD'];
Print;
.store
save Born`typeIU'`typeID'`typeFU'`typeFD'.sav [born`typeIU'`typeID'`typeFU'`typeFD'];
.end
.sort
off stat;
#call SoftBhabha(`typeIU',`typeID',`typeFU',`typeFD',1)
.sort
on stat;
b alpha,pi,volum,prop,Qs,Ts,qel,Lnla,BORN`typeIU'`typeID'`typeFU'`typeFD';
print Soft`typeIU'`typeID'`typeFU'`typeFD';
.end
*---
* form -l -d typeIU=12 -d typeID=12 -d typeFU=12 -d typeFD=12 Brem_Bhabha.frm
*            for e^+ e^- --> e^- e^+
