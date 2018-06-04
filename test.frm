#-
#if 0
 === сокращение дробей ===
sym a,b;
Cfun f;
local A = (a+b)/(a+b);
local B = f(a,b)/f(a,b);
id f?(a?,b?)/f?(a?,b?) = 1;
print;
#endif
#+

* === №1 ===

sym mel,Q,e,volum,SS,UU,TT;
Cfun den;
index mu,nu,smImu;
vector p1,p2,p3,p4;

#call Globals();

local F=       + den(1,mel,Q)^2*e^4*volum * (
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,p1,nu,p4)
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,p1,nu   )*mel
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,p2,nu,p4)
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,p2,nu   )*mel
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,   nu,p4)*mel
          + 1/4*g_(1,nu,p1,smImu,p2,smImu,   nu   )*mel^2
          + 1/4*g_(1,nu,p1,smImu,   smImu,p1,nu,p4)*mel
          + 1/4*g_(1,nu,p1,smImu,   smImu,p1,nu   )*mel^2
          + 1/4*g_(1,nu,p1,smImu,   smImu,p2,nu,p4)*mel
          + 1/4*g_(1,nu,p1,smImu,   smImu,p2,nu   )*mel^2
          + 1/4*g_(1,nu,p1,smImu,   smImu,   nu,p4)*mel^2
          + 1/4*g_(1,nu,p1,smImu,   smImu,   nu   )*mel^3
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,p1,nu,p4)
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,p1,nu   )*mel
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,p2,nu,p4)
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,p2,nu   )*mel
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,   nu,p4)*mel
          + 1/4*g_(1,nu,p2,smImu,p2,smImu,   nu   )*mel^2
          + 1/4*g_(1,nu,p2,smImu,   smImu,p1,nu,p4)*mel
          + 1/4*g_(1,nu,p2,smImu,   smImu,p1,nu   )*mel^2
          + 1/4*g_(1,nu,p2,smImu,   smImu,p2,nu,p4)*mel
          + 1/4*g_(1,nu,p2,smImu,   smImu,p2,nu   )*mel^2
          + 1/4*g_(1,nu,p2,smImu,   smImu,   nu,p4)*mel^2
          + 1/4*g_(1,nu,p2,smImu,   smImu,   nu   )*mel^3
          + 1/4*g_(1,nu,   smImu,p2,smImu,p1,nu,p4)*mel
          + 1/4*g_(1,nu,   smImu,p2,smImu,p1,nu   )*mel^2
          + 1/4*g_(1,nu,   smImu,p2,smImu,p2,nu,p4)*mel
          + 1/4*g_(1,nu,   smImu,p2,smImu,p2,nu   )*mel^2
          + 1/4*g_(1,nu,   smImu,p2,smImu,   nu,p4)*mel^2
          + 1/4*g_(1,nu,   smImu,p2,smImu,   nu   )*mel^3
          + 1/4*g_(1,nu,   smImu,   smImu,p1,nu,p4)*mel^2
          + 1/4*g_(1,nu,   smImu,   smImu,p1,nu   )*mel^3
          + 1/4*g_(1,nu,   smImu,   smImu,p2,nu,p4)*mel^2
          + 1/4*g_(1,nu,   smImu,   smImu,p2,nu   )*mel^3
          + 1/4*g_(1,nu,   smImu,   smImu,   nu,p4)*mel^3
          + 1/4*g_(1,nu,   smImu,   smImu,   nu   )*mel^4
          )

;
trace4 1;
id p1.p1 = 0;
id p2.p2 = mel^2;
id p3.p3 = mel^2;
id p4.p4 = 0;

#define iu  "1"
#define id  "12"
#define fu  "12"
#define fd  "1"

id p1.p2=1/2*(SS-pm(`iu')^2-pm(`id')^2);
id p3.p4=1/2*(SS-pm(`fd')^2-pm(`fu')^2);
id p2.p3=-1/2*(UU-pm(`id')^2-pm(`fd')^2);
id p1.p4=-1/2*(UU-pm(`iu')^2-pm(`fu')^2);
id p2.p4=-1/2*(TT-pm(`id')^2-pm(`fu')^2);
id p1.p3=-1/2*(TT-pm(`iu')^2-pm(`fd')^2);
id TT = pm(`iu')^2+pm(`id')^2+pm(`fu')^2+pm(`fd')^2-SS-UU;
id mgm = 0;

#if 0
id den(1,mel,Q)=1/2/p1.p2;
id SS = mel^2+2*p1.p2;
id UU = mel^2-2*p2.p3;

id p1.p2=1/2*(SS-pm(`iu')^2-pm(`id')^2);
id p3.p4=1/2*(SS-pm(`fd')^2-pm(`fu')^2);
id p2.p3=-1/2*(UU-pm(`id')^2-pm(`fd')^2);
id p1.p4=-1/2*(UU-pm(`iu')^2-pm(`fu')^2);
id p2.p4=-1/2*(TT-pm(`id')^2-pm(`fu')^2);
id p1.p3=-1/2*(TT-pm(`iu')^2-pm(`fd')^2);
id TT = pm(`iu')^2+pm(`id')^2+pm(`fu')^2+pm(`fd')^2-SS-UU;
id mgm = 0;
#endif
    bracket volum,e,mel,den;
    print +s;
    .end :kinematics applied;
