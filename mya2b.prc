#procedure mya2b(a,b)
*------------------
*  no .sort version
*  first generation
*------------------

.sort :mya2bStart;
#include Feynman.h;

* x3,+argument:
*37    id `qup'= `qdn'+1;
*38    id `qdn'= `qup'-1;
*
*35    id `g'= `e'/stw;
*36    id `e'= `g'*stw;
*1     id `mz'^n?!{,0} = (`mw'/ctw)^n;
*2     id `mw'^n?!{,0}= (ctw*`mz')^n;
*3     id `stw'^2= 1 - `ctw'^2;
*4     id `ctw'*`ctw' = 1 - `stw'^2;
*
*5     id gd6(ii?)= 1 + gd5(ii);
*6     id gd5(ii?)= - 1 + gd6(ii);
*7     id gd7(ii?)= 1 - gd5(ii);
*8     id gd5(ii?)= 1 - gd7(ii);
*9     id gd6(ii?)= 2 - gd7(ii);
*10    id gd7(ii?)= 2 - gd6(ii);
*
* x3:
*11    id `vmaup'= `vpaup' - 2*i3up;
*12    id `vpaup'= `vmaup' + 2*i3up;
*13    id `vmadn'= `vpadn' - 2*i3dn;
*14    id `vpadn'= `vmadn' + 2*i3dn;
*15    id `vmaup'= `vpadn' - 2*i3dn - 2 + 2*ctw^2;
*16    id `vpadn'= `vmaup' - 2*i3up + 2 - 2*ctw^2;
*17    id `vmaup'= `vmadn'          - 2 + 2*ctw^2;
*18    id `vmadn'= `vmaup'          + 2 - 2*ctw^2;
*19    id `vpaup'= `vpadn' - 4*i3dn - 2 + 2*ctw^2;
*20    id `vpadn'= `vpaup' - 4*i3up + 2 - 2*ctw^2;
*21    id `vpaup'= `vmadn' - 2*i3dn - 2 + 2*ctw^2;
*22    id `vmadn'= `vpaup' - 2*i3up + 2 - 2*ctw^2;
* x3:
*23    id `vmaen'= `vpaen' - 2*i3en;
*24    id `vpaen'= `vmaen' + 2*i3en;
*25    id `vmael'= `vpael' - 2*i3el;
*26    id `vpael'= `vmael' + 2*i3el;
* x1:
*27    id `vmaen'= `vpael' - 2*i3el - 2 + 2*ctw^2;
*28    id `vpael'= `vmaen' - 2*i3en + 2 - 2*ctw^2;
*29    id `vmaen'= `vmael'          - 2 + 2*ctw^2;
*30    id `vmael'= `vmaen'          + 2 - 2*ctw^2;
*31    id `vpaen'= `vpael' - 4*i3el - 2 + 2*ctw^2;
*32    id `vpael'= `vpaen' - 4*i3en + 2 - 2*ctw^2;
*33    id `vpaen'= `vmael' - 2*i3el - 2 + 2*ctw^2;
*34    id `vmael'= `vpaen' - 2*i3en + 2 - 2*ctw^2;
*
* x3:
*39    id `vmaup'= `vup' - i3up;
*40    id `vpaup'= `vup' + i3up;
*41    id `vmadn'= `vdn' - i3dn;
*42    id `vpadn'= `vdn' + i3dn;
*43    id `vmael'= `vel' - i3el;
*44    id `vpael'= `vel' + i3el;
*45 unused
* x1:
*46    id `vpaen'= 'i3en' + i3en;

*------------------ 1-10
#$s1L = `a'-mz;
#$s1R = `b'-mw;
#if ("`$s1L'"==0) && ("`$s1R'"==0)
    id `a'^n?!{,0} = (`b'/ctw)^n;
#endif

#$s2L = `a'-mw;
#$s2R = `b'-mz;
#if ("`$s2L'"==0) && ("`$s2R'"==0)
    id `a'^n?!{,0}= (ctw*`b')^n;
#endif


#$s3L = `a'-stw;
#$s3R = `b'-ctw;
#if ("`$s3L'"==0) && ("`$s3R'"==0)
    id `a'^2= 1 - `b'^2;
    id stw^-2*ctw^-1= ctw/stw^2 + 1/ctw;
#endif

#$s4L = `a'-ctw;
#$s4R = `b'-stw;
#if ("`$s4L'"==0) && ("`$s4R'"==0)
    id `a'*`a' = 1 - `b'^2;
    id ctw^-2*stw^-1= stw/ctw^2 + 1/stw;
    id ctw^-4=(1+stw^2/ctw^2)/ctw^2;
    id ctw^-2=1+stw^2/ctw^2;
#endif


#$s5L = `a'-gd6;
#$s5R = `b'-gd5;
#if ("`$s5L'"==0) && ("`$s5R'"==0)
    id gd6(ii?)= 1 + gd5(ii);
#endif

#$s6L = `a'-gd5;
#$s6R = `b'-gd6;
#if ("`$s6L'"==0) && ("`$s6R'"==0)
    id gd5(ii?)= - 1 + gd6(ii);
#endif


#$s7L = `a'-gd7;
#$s7R = `b'-gd5;
#if ("`$s7L'"==0) && ("`$s7R'"==0)
    id gd7(ii?)= 1 - gd5(ii);
#endif

#$s8L = `a'-gd5;
#$s8R = `b'-gd7;
#if ("`$s8L'"==0) && ("`$s8R'"==0)
    id gd5(ii?)= 1 - gd7(ii);
#endif


#$s9L = `a'-gd6;
#$s9R = `b'-gd7;
#if ("`$s9L'"==0) && ("`$s9R'"==0)
    id gd6(ii?)= 2 - gd7(ii);
#endif

#$s10L= `a'-gd7;
#$s10R= `b'-gd6;
#if ("`$s10L'"==0) && ("`$s10R'"==0)
    id gd7(ii?)= 2 - gd6(ii);
#endif

*--------------------------------
*--------------------------------11
#$s11L= `a'-vmaup;
#$s11R= `b'-vpaup;
#if ("`$s11L'"==0)  && ("`$s11R'"==0)
    id `a' = `b' - 2*i3up;
#endif;
#$s211L= `a'-vmach;
#$s211R= `b'-vpach;
#if ("`$s211L'"==0) && ("`$s211R'"==0)
    id `a' = `b' - 2*i3ch;
#endif;
#$s311L= `a'-vmatp;
#$s311R= `b'-vpatp;
#if ("`$s311L'"==0) && ("`$s311R'"==0)
    id `a' = `b' - 2*i3tp;
#endif;
*--------------------------------12
#$s12L= `a'-vpaup;
#$s12R= `b'-vmaup;
#if ("`$s12L'"==0)  && ("`$s12R'"==0)
    id `a'= `b' + 2*i3up;
#endif
#$s212L= `a'-vpach;
#$s212R= `b'-vmach;
#if ("`$s212L'"==0) && ("`$s212R'"==0)
    id `a'= `b' + 2*i3ch;
#endif
#$s312L= `a'-vpatp;
#$s312R= `b'-vmatp;
#if ("`$s312L'"==0) && ("`$s312R'"==0)
    id `a'= `b' + 2*i3tp;
#endif
*--------------------------------13
#$s13L= `a'-vmadn;
#$s13R= `b'-vpadn;
#if ("`$s13L'"==0)  && ("`$s13R'"==0)
    id `a'= `b' - 2*i3dn;
#endif
#$s213L= `a'-vmast;
#$s213R= `b'-vpast;
#if ("`$s213L'"==0) && ("`$s213R'"==0)
    id `a'= `b' - 2*i3st;
#endif
#$s313L= `a'-vmabt;
#$s313R= `b'-vpabt;
#if ("`$s313L'"==0) && ("`$s313R'"==0)
    id `a'= `b' - 2*i3bt;
#endif
*--------------------------------14
#$s14L= `a'-vpadn;
#$s14R= `b'-vmadn;
#if ("`$s14L'"==0)  && ("`$s14R'"==0)
    id `a'= `b' + 2*i3dn;
#endif
#$s214L= `a'-vpast;
#$s214R= `b'-vmast;
#if ("`$s214L'"==0) && ("`$s214R'"==0)
    id `a'= `b' + 2*i3st;
#endif
#$s314L= `a'-vpabt;
#$s314R= `b'-vmabt;
#if ("`$s314L'"==0) && ("`$s314R'"==0)
    id `a'= `b' + 2*i3bt;
#endif
*--------------------------------15
#$s15L= `a'-vmaup;
#$s15R= `b'-vpadn;
#if ("`$s15L'"==0)  && ("`$s15R'"==0)
    id `a'= `b' - 2*i3dn - 2 + 2*ctw^2;
#endif
#$s215L= `a'-vmach;
#$s215R= `b'-vpast;
#if ("`$s215L'"==0) && ("`$s215R'"==0)
    id `a'= `b' - 2*i3st - 2 + 2*ctw^2;
#endif
#$s315L= `a'-vmatp;
#$s315R= `b'-vpabt;
#if ("`$s315L'"==0) && ("`$s315R'"==0)
    id `a'= `b' - 2*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------16
#$s16L= `a'-vpadn;
#$s16R= `b'-vmaup;
#if ("`$s16L'"==0)  && ("`$s16R'"==0)
    id `a'= `b' - 2*i3up + 2 - 2*ctw^2;
#endif
#$s216L= `a'-vpast;
#$s216R= `b'-vmach;
#if ("`$s216L'"==0) && ("`$s216R'"==0)
    id `a'= `b' - 2*i3ch + 2 - 2*ctw^2;
#endif
#$s316L= `a'-vpabt;
#$s316R= `b'-vmatp;
#if ("`$s316L'"==0) && ("`$s316R'"==0)
    id `a'= `b' - 2*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------17
#$s17L= `a'-vmaup;
#$s17R= `b'-vmadn;
#if ("`$s17L'"==0)  && ("`$s17R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
#$s217L= `a'-vmach;
#$s217R= `b'-vmast;
#if ("`$s217L'"==0) && ("`$s217R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
#$s317L= `a'-vmatp;
#$s317R= `b'-vmabt;
#if ("`$s317L'"==0) && ("`$s317R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
*--------------------------------18
#$s18L= `a'-vmadn;
#$s18R= `b'-vmaup;
#if ("`$s18L'"==0)  && ("`$s18R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
#$s218L= `a'-vmast;
#$s218R= `b'-vmach;
#if ("`$s218L'"==0) && ("`$s218R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
#$s318L= `a'-vmabt;
#$s318R= `b'-vmatp;
#if ("`$s318L'"==0) && ("`$s318R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
*--------------------------------19
#$s19L= `a'-vpaup;
#$s19R= `b'-vpadn;
#if ("`$s19L'"==0)  && ("`$s19R'"==0)
   id `a'= `b' - 4*i3dn - 2 + 2*ctw^2;
#endif
#$s219L= `a'-vpach;
#$s219R= `b'-vpast;
#if ("`$s219L'"==0) && ("`$s219R'"==0)
   id `a'= `b' - 4*i3st - 2 + 2*ctw^2;
#endif
#$s319L= `a'-vpatp;
#$s319R= `b'-vpabt;
#if ("`$s319L'"==0) && ("`$s319R'"==0)
   id `a'= `b' - 4*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------20
#$s20L= `a'-vpadn;
#$s20R= `b'-vpaup;
#if ("`$s20L'"==0)  && ("`$s20R'"==0)
    id `a'= `b' - 4*i3up + 2 - 2*ctw^2;
#endif
#$s220L= `a'-vpast;
#$s220R= `b'-vpach;
#if ("`$s220L'"==0) && ("`$s220R'"==0)
    id `a'= `b' - 4*i3ch + 2 - 2*ctw^2;
#endif
#$s320L= `a'-vpabt;
#$s320R= `b'-vpatp;
#if ("`$s320L'"==0) && ("`$s320R'"==0)
    id `a'= `b' - 4*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------21
#$s21L= `a'-vpaup;
#$s21R= `b'-vmadn;
#if ("`$s21L'"==0)  && ("`$s21R'"==0)
    id `a'= `b' - 2*i3dn - 2 + 2*ctw^2;
#endif
#$s221L= `a'-vpach;
#$s221R= `b'-vmast;
#if ("`$s221L'"==0) && ("`$s221R'"==0)
    id `a'= `b' - 2*i3st - 2 + 2*ctw^2;
#endif
#$s321L= `a'-vpatp;
#$s321R= `b'-vmabt;
#if ("`$s321L'"==0) && ("`$s321R'"==0)
    id `a'= `b' - 2*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------22
#$s22L= `a'-vmadn;
#$s22R= `b'-vpaup;
#if ("`$s22L'"==0)  && ("`$s22R'"==0)
    id `a'= `b' - 2*i3up + 2 - 2*ctw^2;
#endif
#$s222L= `a'-vmast;
#$s222R= `b'-vpach;
#if ("`$s222L'"==0) && ("`$s222R'"==0)
    id `a'= `b' - 2*i3ch + 2 - 2*ctw^2;
#endif
#$s322L= `a'-vmabt;
#$s322R= `b'-vpatp;
#if ("`$s322L'"==0) && ("`$s322R'"==0)
    id `a'= `b' - 2*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------23
#$s23L= `a'-vmaen;
#$s23R= `b'-vpaen;
#if ("`$s23L'"==0)  && ("`$s23R'"==0)
    id `a' = `b' - 2*i3en;
#endif;
#$s223L= `a'-vmamn;
#$s223R= `b'-vpamn;
#if ("`$s223L'"==0) && ("`$s223R'"==0)
    id `a' = `b' - 2*i3mn;
#endif;
#$s323L= `a'-vmatn;
#$s323R= `b'-vpatn;
#if ("`$s323L'"==0) && ("`$s323R'"==0)
    id `a' = `b' - 2*i3tn;
#endif;
*--------------------------------24
#$s24L= `a'-vpaen;
#$s24R= `b'-vmaen;
#if ("`$s24L'"==0)  && ("`$s24R'"==0)
    id `a'= `b' + 2*i3en;
#endif
#$s224L= `a'-vpamn;
#$s224R= `b'-vmamn;
#if ("`$s224L'"==0) && ("`$s224R'"==0)
    id `a'= `b' + 2*i3mn;
#endif
#$s324L= `a'-vpatn;
#$s324R= `b'-vmatn;
#if ("`$s324L'"==0) && ("`$s324R'"==0)
    id `a'= `b' + 2*i3tn;
#endif
*--------------------------------25
#$s25L= `a'-vmael;
#$s25R= `b'-vpael;
#if ("`$s25L'"==0)  && ("`$s25R'"==0)
    id `a'= `b' - 2*i3el;
#endif
#$s225L= `a'-vmamo;
#$s225R= `b'-vpamo;
#if ("`$s225L'"==0) && ("`$s225R'"==0)
    id `a'= `b' - 2*i3mo;
#endif
#$s325L= `a'-vmata;
#$s325R= `b'-vpata;
#if ("`$s325L'"==0) && ("`$s325R'"==0)
    id `a'= `b' - 2*i3ta;
#endif
*--------------------------------26
#$s26L= `a'-vpael;
#$s26R= `b'-vmael;
#if ("`$s26L'"==0)  && ("`$s26R'"==0)
    id `a'= `b' + 2*i3el;
#endif
#$s226L= `a'-vpamo;
#$s226R= `b'-vmamo;
#if ("`$s226L'"==0) && ("`$s226R'"==0)
    id `a'= `b' + 2*i3mo;
#endif
#$s326L= `a'-vpata;
#$s326R= `b'-vmata;
#if ("`$s326L'"==0) && ("`$s326R'"==0)
    id `a'= `b' + 2*i3ta;
#endif
*--------------------------------
*--------------------------------27-36

#$s27L= `a'-vmaen;
#$s27R= `b'-vpael;
#if ("`$s27L'"==0) && ("`$s27R'"==0)
    id `a'= `b' - 2*i3el - 2 + 2*ctw^2;
#endif

#$s28L= `a'-vpael;
#$s28R= `b'-vmaen;
#if ("`$s28L'"==0) && ("`$s28R'"==0)
    id `a'= `b' - 2*i3en + 2 - 2*ctw^2;
#endif
    
#$s29L= `a'-vmaen;
#$s29R= `b'-vmael;
#if ("`$s29L'"==0) && ("`$s29R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif

#$s30L= `a'-vmael;
#$s30R= `b'-vmaen;
#if ("`$s30L'"==0) && ("`$s30R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif

#$s31L= `a'-vpaen;
#$s31R= `b'-vpael;
#if ("`$s31L'"==0) && ("`$s31R'"==0)
   id `a'= `b' - 4*i3el - 2 + 2*ctw^2;
#endif

#$s32L= `a'-vpael;
#$s32R= `b'-vpaen;
#if ("`$s32L'"==0) && ("`$s32R'"==0)
    id `a'= `b' - 4*i3en + 2 - 2*ctw^2;
#endif

#$s33L= `a'-vpaen;
#$s33R= `b'-vmael;
#if ("`$s33L'"==0) && ("`$s33R'"==0)
    id `a'= `b' - 2*i3el - 2 + 2*ctw^2;
#endif

#$s34L= `a'-vmael;
#$s34R= `b'-vpaen;
#if ("`$s34L'"==0) && ("`$s34R'"==0)
    id `a'= `b' - 2*i3en + 2 - 2*ctw^2;
#endif

#$s35L= `a'-g;
#$s35R= `b'-e;
#if ("`$s35L'"==0) && ("`$s35R'"==0)
    id `a'= `b'/stw;
#endif
#$s36L= `a'-e;
#$s36R= `b'-g;
#if ("`$s36L'"==0) && ("`$s36R'"==0)
    id `a'= `b'*stw;
#endif

*--------------------------------
*--------------------------------37
#$s37L= `a'-qup;
#$s37R= `b'-qdn;
#if ("`$s37L'"==0)  && ("`$s37R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
#$s237L= `a'-qch;
#$s237R= `b'-qst;
#if ("`$s237L'"==0) && ("`$s237R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
#$s337L= `a'-qtp;
#$s337R= `b'-qbt;
#if ("`$s337L'"==0) && ("`$s337R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
*--------------------------------38
#$s38L= `a'-qdn;
#$s38R= `b'-qup;
#if ("`$s38L'"==0)  && ("`$s38R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
#$s238L= `a'-qst;
#$s238R= `b'-qch;
#if ("`$s238L'"==0) && ("`$s238R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
#$s338L= `a'-qbt;
#$s338R= `b'-qtp;
#if ("`$s338L'"==0) && ("`$s338R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
*--------------------------------39
#$s39L= `a'-vmaup;
#$s39R= `b'-vup;
#if ("`$s39L'"==0)  && ("`$s39R'"==0)
    id `a'= `b' - i3up;
#endif
#$s239L= `a'-vmach;
#$s239R= `b'-vch;
#if ("`$s239L'"==0) && ("`$s239R'"==0)
    id `a'= `b' - i3ch;
#endif
#$s339L= `a'-vmatp;
#$s339R= `b'-vtp;
#if ("`$s339L'"==0) && ("`$s339R'"==0)
    id `a'= `b' - i3tp;
#endif
*--------------------------------40
#$s40L= `a'-vpaup;
#$s40R= `b'-vup;
#if ("`$s40L'"==0)  && ("`$s40R'"==0)
    id `a'= `b' + i3up;
#endif
#$s240L= `a'-vpach;
#$s240R= `b'-vch;
#if ("`$s240L'"==0) && ("`$s240R'"==0)
    id `a'= `b' + i3ch;
#endif
#$s340L= `a'-vpatp;
#$s340R= `b'-vtp;
#if ("`$s340L'"==0) && ("`$s340R'"==0)
    id `a'= `b' + i3tp;
#endif
*--------------------------------41
#$s41L= `a'-vmadn;
#$s41R= `b'-vdn;
#if ("`$s41L'"==0)  && ("`$s41R'"==0)
    id `a'= `b' - i3dn;
#endif
#$s241L= `a'-vmast;
#$s241R= `b'-vst;
#if ("`$s241L'"==0) && ("`$s241R'"==0)
    id `a'= `b' - i3st;
#endif
#$s341L= `a'-vmabt;
#$s341R= `b'-vbt;
#if ("`$s341L'"==0) && ("`$s341R'"==0)
    id `a'= `b' - i3bt;
#endif
*--------------------------------42
#$s42L= `a'-vpadn;
#$s42R= `b'-vdn;
#if ("`$s42L'"==0)  && ("`$s42R'"==0)
    id `a'= `b' + i3dn;
#endif
#$s242L= `a'-vpast;
#$s242R= `b'-vst;
#if ("`$s242L'"==0) && ("`$s242R'"==0)
    id `a'= `b' + i3st;
#endif
#$s342L= `a'-vpabt;
#$s342R= `b'-vbt;
#if ("`$s342L'"==0) && ("`$s342R'"==0)
    id `a'= `b' + i3bt;
#endif
*--------------------------------43
#$s43L= `a'-vmael;
#$s43R= `b'-vel;
#if ("`$s43L'"==0)  && ("`$s43R'"==0)
    id `a'= `b' - i3el;
#endif
#$s243L= `a'-vmamo;
#$s243R= `b'-vmo;
#if ("`$s243L'"==0) && ("`$s243R'"==0)
    id `a'= `b' - i3mo;
#endif
#$s343L= `a'-vmata;
#$s343R= `b'-vta;
#if ("`$s343L'"==0) && ("`$s343R'"==0)
    id `a'= `b' - i3ta;
#endif
*--------------------------------44
#$s44L= `a'-vpael;
#$s44R= `b'-vel;
#if ("`$s44L'"==0)  && ("`$s44R'"==0)
    id `a'= `b' + i3el;
#endif
#$s244L= `a'-vpamo;
#$s244R= `b'-vmo;
#if ("`$s244L'"==0) && ("`$s244R'"==0)
    id `a'= `b' + i3mo;
#endif
#$s344L= `a'-vpata;
#$s344R= `b'-vta;
#if ("`$s344L'"==0) && ("`$s344R'"==0)
    id `a'= `b' + i3ta;
#endif
*--------------------------------
*45 unused
*--------------------------------46

#$s46L= `a'-vpaen;
#$s46R= `b'-i3en;
#if ("`$s46L'"==0) && ("`$s46R'"==0)
    id `a'= 'b' + i3en;
#endif
#endprocedure
*------------
