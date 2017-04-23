#procedure BornBhabha(iu,id,fu,fd)
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
.sort
g Born`iu'`id'`fu'`fd' =
#do vb={1,2}
   + i_*Vb(ii,p1,h1)*vert(-`vb',`id',-`iu',nu,ii)*U(ii,p2,h2)*pr(`vb',nu,mu,Q)
       *Ub(jj,p3,h3)*vert(`vb',`fu',-`fd',mu,jj)*V(jj,p4,h4)
   - i_*Ub(ii2,p3,h3)*vert(-`vb',`id',-`fd',mu1,ii2)*U(ii2,p2,h2)*pr(`vb',nu1,mu1,P)
       *Vb(ii1,p1,h1)*vert(`vb',`fu',-`iu',nu1,ii1)*V(ii1,p4,h4)
#enddo
;
print +s;
.sort:born created;

#call FeynmanRules(1)
print +s;
.sort:FeynmanRules applied;

id gd(ii,Q) = gd(ii,p1)+gd(ii,p2);
id gd(jj,Q) = gd(jj,p3)+gd(jj,p4);
id gd(ii1,P)=-gd(ii1,p1)+gd(ii1,p4);
id gd(ii2,P)= gd(ii2,p2)-gd(ii2,p3);
print +s;
.sort:Q and P expanded;

#call DiracEquation(Born`iu'`id'`fu'`fd',-1)
#call Convert(1)
#call a2b(gd7,gd6)
id den(1,mh,Q)=prop(mh^2,Qs);
id den(1,mh,P)=prop(mh^2,Ts);
id g^2*stw^2*den(1,0,Q)=4*pi*alpha/Qs;
id g^2*stw^2*den(1,0,P)=4*pi*alpha/Ts;
id g^2/ctw^2*den(1,mp?,Q)=4*pi*alpha/Qs*4*chi(mp^2,s);
id g^2/ctw^2*den(1,mp?,P)=4*pi*alpha/Ts*4*chi(mp^2,t);
.sort
#call MakeAmpSquare(Born`iu'`id'`fu'`fd',amplitudeSquared,1/4*volum)
#call Trace(amplitudeSquared,[born`iu'`id'`fu'`fd'],1,1)
.sort
drop amplitudeSquared;
.sort
*id volum=1/32/pi/s*d(cos);    
id p1.p2=1/2*(Qs+pm(`iu')^2+pm(`id')^2);
id p3.p4=1/2*(Qs+pm(`fd')^2+pm(`fu')^2);
id p2.p3=-1/2*(Ts+pm(`id')^2+pm(`fd')^2);
id p1.p4=-1/2*(Ts+pm(`iu')^2+pm(`fu')^2);
id p2.p4=-1/2*(Us+pm(`id')^2+pm(`fu')^2);
id p1.p3=-1/2*(Us+pm(`iu')^2+pm(`fd')^2);
.sort
#endprocedure
*------------
