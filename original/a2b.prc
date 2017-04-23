#procedure a2b(a,b)
*------------------
*  no .sort version
*  first generation
*------------------
#$s1L = `a'-mz;
#$s1R = `b'-mw;
#$s2L = `a'-mw;
#$s2R = `b'-mz;

#$s3L = `a'-stw;
#$s3R = `b'-ctw;
#$s4L = `a'-ctw;
#$s4R = `b'-stw;

#$s5L = `a'-gd6;
#$s5R = `b'-gd5;
#$s6L = `a'-gd5;
#$s6R = `b'-gd6;

#$s7L = `a'-gd7;
#$s7R = `b'-gd5;
#$s8L = `a'-gd5;
#$s8R = `b'-gd7;

#$s9L = `a'-gd6;
#$s9R = `b'-gd7;
#$s10L= `a'-gd7;
#$s10R= `b'-gd6;

#$s11L= `a'-vmaup;
#$s11R= `b'-vpaup;
#$s12L= `a'-vpaup;
#$s12R= `b'-vmaup;

#$s13L= `a'-vmadn;
#$s13R= `b'-vpadn;
#$s14L= `a'-vpadn;
#$s14R= `b'-vmadn;

#$s15L= `a'-vmaup;
#$s15R= `b'-vpadn;
#$s16L= `a'-vpadn;
#$s16R= `b'-vmaup;

#$s17L= `a'-vmaup;
#$s17R= `b'-vmadn;
#$s18L= `a'-vmadn;
#$s18R= `b'-vmaup;

#$s19L= `a'-vpaup;
#$s19R= `b'-vpadn;
#$s20L= `a'-vpadn;
#$s20R= `b'-vpaup;

#$s21L= `a'-vpaup;
#$s21R= `b'-vmadn;
#$s22L= `a'-vmadn;
#$s22R= `b'-vpaup;

#$s23L= `a'-vmaen;
#$s23R= `b'-vpaen;
#$s24L= `a'-vpaen;
#$s24R= `b'-vmaen;

#$s25L= `a'-vmael;
#$s25R= `b'-vpael;
#$s26L= `a'-vpael;
#$s26R= `b'-vmael;

#$s27L= `a'-vmaen;
#$s27R= `b'-vpael;
#$s28L= `a'-vpael;
#$s28R= `b'-vmaen;

#$s29L= `a'-vmaen;
#$s29R= `b'-vmael;
#$s30L= `a'-vmael;
#$s30R= `b'-vmaen;

#$s31L= `a'-vpaen;
#$s31R= `b'-vpael;
#$s32L= `a'-vpael;
#$s32R= `b'-vpaen;

#$s33L= `a'-vpaen;
#$s33R= `b'-vmael;
#$s34L= `a'-vmael;
#$s34R= `b'-vpaen;

#$s35L= `a'-g;
#$s35R= `b'-e;
#$s36L= `a'-e;
#$s36R= `b'-g;

#$s37L= `a'-qup;
#$s37R= `b'-qdn;
#$s38L= `a'-qdn;
#$s38R= `b'-qup;

#$s39L= `a'-vmaup;
#$s39R= `b'-vup;
#$s40L= `a'-vpaup;
#$s40R= `b'-vup;

#$s41L= `a'-vmadn;
#$s41R= `b'-vdn;
#$s42L= `a'-vpadn;
#$s42R= `b'-vdn;

#$s43L= `a'-vmael;
#$s43R= `b'-vel;
#$s44L= `a'-vpael;
#$s44R= `b'-vel;

*45 unused

#$s46L= `a'-vpaen;
#$s46R= `b'-i3en;
*
*------------------
* second generation
#$s211L= `a'-vmach;
#$s211R= `b'-vpach;
#$s212L= `a'-vpach;
#$s212R= `b'-vmach;
#$s213L= `a'-vmast;
#$s213R= `b'-vpast;
#$s214L= `a'-vpast;
#$s214R= `b'-vmast;
#$s215L= `a'-vmach;
#$s215R= `b'-vpast;
#$s216L= `a'-vpast;
#$s216R= `b'-vmach;
#$s217L= `a'-vmach;
#$s217R= `b'-vmast;
#$s218L= `a'-vmast;
#$s218R= `b'-vmach;
#$s219L= `a'-vpach;
#$s219R= `b'-vpast;
#$s220L= `a'-vpast;
#$s220R= `b'-vpach;
#$s221L= `a'-vpach;
#$s221R= `b'-vmast;
#$s222L= `a'-vmast;
#$s222R= `b'-vpach;
#$s223L= `a'-vmamn;
#$s223R= `b'-vpamn;
#$s224L= `a'-vpamn;
#$s224R= `b'-vmamn;
#$s225L= `a'-vmamo;
#$s225R= `b'-vpamo;
#$s226L= `a'-vpamo;
#$s226R= `b'-vmamo;
*
#$s237L= `a'-qch;
#$s237R= `b'-qst;
#$s238L= `a'-qst;
#$s238R= `b'-qch;
#$s239L= `a'-vmach;
#$s239R= `b'-vch;
#$s240L= `a'-vpach;
#$s240R= `b'-vch;
#$s241L= `a'-vmast;
#$s241R= `b'-vst;
#$s242L= `a'-vpast;
#$s242R= `b'-vst;
#$s243L= `a'-vmamo;
#$s243R= `b'-vmo;
#$s244L= `a'-vpamo;
#$s244R= `b'-vmo;
*------------------
* third  generation
#$s311L= `a'-vmatp;
#$s311R= `b'-vpatp;
#$s312L= `a'-vpatp;
#$s312R= `b'-vmatp;
#$s313L= `a'-vmabt;
#$s313R= `b'-vpabt;
#$s314L= `a'-vpabt;
#$s314R= `b'-vmabt;
#$s315L= `a'-vmatp;
#$s315R= `b'-vpabt;
#$s316L= `a'-vpabt;
#$s316R= `b'-vmatp;
#$s317L= `a'-vmatp;
#$s317R= `b'-vmabt;
#$s318L= `a'-vmabt;
#$s318R= `b'-vmatp;
#$s319L= `a'-vpatp;
#$s319R= `b'-vpabt;
#$s320L= `a'-vpabt;
#$s320R= `b'-vpatp;
#$s321L= `a'-vpatp;
#$s321R= `b'-vmabt;
#$s322L= `a'-vmabt;
#$s322R= `b'-vpatp;
#$s323L= `a'-vmatn;
#$s323R= `b'-vpatn;
#$s324L= `a'-vpatn;
#$s324R= `b'-vmatn;
#$s325L= `a'-vmata;
#$s325R= `b'-vpata;
#$s326L= `a'-vpata;
#$s326R= `b'-vmata;
*
#$s337L= `a'-qtp;
#$s337R= `b'-qbt;
#$s338L= `a'-qbt;
#$s338R= `b'-qtp;
#$s339L= `a'-vmatp;
#$s339R= `b'-vtp;
#$s340L= `a'-vpatp;
#$s340R= `b'-vtp;
#$s341L= `a'-vmabt;
#$s341R= `b'-vbt;
#$s342L= `a'-vpabt;
#$s342R= `b'-vbt;
#$s343L= `a'-vmata;
#$s343R= `b'-vta;
#$s344L= `a'-vpata;
#$s344R= `b'-vta;
*------------------
*
#if ("`$s1L'"==0) && ("`$s1R'"==0)
    id `a'^n?!{,0} = (`b'/ctw)^n;
#endif

#if ("`$s2L'"==0) && ("`$s2R'"==0)
    id `a'^n?!{,0}= (ctw*`b')^n;
#endif

#if ("`$s3L'"==0) && ("`$s3R'"==0)
    id `a'^2= 1 - `b'^2;
    id stw^-2*ctw^-1= ctw/stw^2 + 1/ctw;
#endif

#if ("`$s4L'"==0) && ("`$s4R'"==0)
    id `a'*`a' = 1 - `b'^2;
    id ctw^-2*stw^-1= stw/ctw^2 + 1/stw;
    id ctw^-4=(1+stw^2/ctw^2)/ctw^2;
    id ctw^-2=1+stw^2/ctw^2;
#endif

#if ("`$s5L'"==0) && ("`$s5R'"==0)
    id gd6(ii?)= 1 + gd5(ii);
#endif

#if ("`$s6L'"==0) && ("`$s6R'"==0)
    id gd5(ii?)= - 1 + gd6(ii);
#endif

#if ("`$s7L'"==0) && ("`$s7R'"==0)
    id gd7(ii?)= 1 - gd5(ii);
#endif

#if ("`$s8L'"==0) && ("`$s8R'"==0)
    id gd5(ii?)= 1 - gd7(ii);
#endif

#if ("`$s9L'"==0) && ("`$s9R'"==0)
    id gd6(ii?)= 2 - gd7(ii);
#endif

#if ("`$s10L'"==0) && ("`$s10R'"==0)
    id gd7(ii?)= 2 - gd6(ii);
#endif
*--------------------------------
#if ("`$s11L'"==0)  && ("`$s11R'"==0)
    id `a' = `b' - 2*i3up;
#endif;
#if ("`$s211L'"==0) && ("`$s211R'"==0)
    id `a' = `b' - 2*i3ch;
#endif;
#if ("`$s311L'"==0) && ("`$s311R'"==0)
    id `a' = `b' - 2*i3tp;
#endif;
*--------------------------------
#if ("`$s12L'"==0)  && ("`$s12R'"==0)
    id `a'= `b' + 2*i3up;
#endif
#if ("`$s212L'"==0) && ("`$s212R'"==0)
    id `a'= `b' + 2*i3ch;
#endif
#if ("`$s312L'"==0) && ("`$s312R'"==0)
    id `a'= `b' + 2*i3tp;
#endif
*--------------------------------
#if ("`$s13L'"==0)  && ("`$s13R'"==0)
    id `a'= `b' - 2*i3dn;
#endif
#if ("`$s213L'"==0) && ("`$s213R'"==0)
    id `a'= `b' - 2*i3st;
#endif
#if ("`$s313L'"==0) && ("`$s313R'"==0)
    id `a'= `b' - 2*i3bt;
#endif
*--------------------------------
#if ("`$s14L'"==0)  && ("`$s14R'"==0)
    id `a'= `b' + 2*i3dn;
#endif
#if ("`$s214L'"==0) && ("`$s214R'"==0)
    id `a'= `b' + 2*i3st;
#endif
#if ("`$s314L'"==0) && ("`$s314R'"==0)
    id `a'= `b' + 2*i3bt;
#endif
*--------------------------------
#if ("`$s15L'"==0)  && ("`$s15R'"==0)
    id `a'= `b' - 2*i3dn - 2 + 2*ctw^2;
#endif
#if ("`$s215L'"==0) && ("`$s215R'"==0)
    id `a'= `b' - 2*i3st - 2 + 2*ctw^2;
#endif
#if ("`$s315L'"==0) && ("`$s315R'"==0)
    id `a'= `b' - 2*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------
#if ("`$s16L'"==0)  && ("`$s16R'"==0)
    id `a'= `b' - 2*i3up + 2 - 2*ctw^2;
#endif
#if ("`$s216L'"==0) && ("`$s216R'"==0)
    id `a'= `b' - 2*i3ch + 2 - 2*ctw^2;
#endif
#if ("`$s316L'"==0) && ("`$s316R'"==0)
    id `a'= `b' - 2*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------
#if ("`$s17L'"==0)  && ("`$s17R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
#if ("`$s217L'"==0) && ("`$s217R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
#if ("`$s317L'"==0) && ("`$s317R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif
*--------------------------------
#if ("`$s18L'"==0)  && ("`$s18R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
#if ("`$s218L'"==0) && ("`$s218R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
#if ("`$s318L'"==0) && ("`$s318R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif
*--------------------------------
#if ("`$s19L'"==0)  && ("`$s19R'"==0)
   id `a'= `b' - 4*i3dn - 2 + 2*ctw^2;
#endif
#if ("`$s219L'"==0) && ("`$s219R'"==0)
   id `a'= `b' - 4*i3st - 2 + 2*ctw^2;
#endif
#if ("`$s319L'"==0) && ("`$s319R'"==0)
   id `a'= `b' - 4*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------
#if ("`$s20L'"==0)  && ("`$s20R'"==0)
    id `a'= `b' - 4*i3up + 2 - 2*ctw^2;
#endif
#if ("`$s220L'"==0) && ("`$s220R'"==0)
    id `a'= `b' - 4*i3ch + 2 - 2*ctw^2;
#endif
#if ("`$s320L'"==0) && ("`$s320R'"==0)
    id `a'= `b' - 4*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------
#if ("`$s21L'"==0)  && ("`$s21R'"==0)
    id `a'= `b' - 2*i3dn - 2 + 2*ctw^2;
#endif
#if ("`$s221L'"==0) && ("`$s221R'"==0)
    id `a'= `b' - 2*i3st - 2 + 2*ctw^2;
#endif
#if ("`$s321L'"==0) && ("`$s321R'"==0)
    id `a'= `b' - 2*i3bt - 2 + 2*ctw^2;
#endif
*--------------------------------
#if ("`$s22L'"==0)  && ("`$s22R'"==0)
    id `a'= `b' - 2*i3up + 2 - 2*ctw^2;
#endif
#if ("`$s222L'"==0) && ("`$s222R'"==0)
    id `a'= `b' - 2*i3ch + 2 - 2*ctw^2;
#endif
#if ("`$s322L'"==0) && ("`$s322R'"==0)
    id `a'= `b' - 2*i3tp + 2 - 2*ctw^2;
#endif
*--------------------------------
#if ("`$s23L'"==0)  && ("`$s23R'"==0)
    id `a' = `b' - 2*i3en;
#endif;
#if ("`$s223L'"==0) && ("`$s223R'"==0)
    id `a' = `b' - 2*i3mn;
#endif;
#if ("`$s323L'"==0) && ("`$s323R'"==0)
    id `a' = `b' - 2*i3tn;
#endif;
*--------------------------------
#if ("`$s24L'"==0)  && ("`$s24R'"==0)
    id `a'= `b' + 2*i3en;
#endif
#if ("`$s224L'"==0) && ("`$s224R'"==0)
    id `a'= `b' + 2*i3mn;
#endif
#if ("`$s324L'"==0) && ("`$s324R'"==0)
    id `a'= `b' + 2*i3tn;
#endif
*--------------------------------
#if ("`$s25L'"==0)  && ("`$s25R'"==0)
    id `a'= `b' - 2*i3el;
#endif
#if ("`$s225L'"==0) && ("`$s225R'"==0)
    id `a'= `b' - 2*i3mo;
#endif
#if ("`$s325L'"==0) && ("`$s325R'"==0)
    id `a'= `b' - 2*i3ta;
#endif
*--------------------------------
#if ("`$s26L'"==0)  && ("`$s26R'"==0)
    id `a'= `b' + 2*i3el;
#endif
#if ("`$s226L'"==0) && ("`$s226R'"==0)
    id `a'= `b' + 2*i3mo;
#endif
#if ("`$s326L'"==0) && ("`$s326R'"==0)
    id `a'= `b' + 2*i3ta;
#endif
*--------------------------------
#if ("`$s27L'"==0) && ("`$s27R'"==0)
    id `a'= `b' - 2*i3el - 2 + 2*ctw^2;
#endif

#if ("`$s28L'"==0) && ("`$s28R'"==0)
    id `a'= `b' - 2*i3en + 2 - 2*ctw^2;
#endif
    
#if ("`$s29L'"==0) && ("`$s29R'"==0)
    id `a'= `b'- 2 + 2*ctw^2;
#endif

#if ("`$s30L'"==0) && ("`$s30R'"==0)
    id `a'= `b'+ 2 - 2*ctw^2;
#endif

#if ("`$s31L'"==0) && ("`$s31R'"==0)
   id `a'= `b' - 4*i3el - 2 + 2*ctw^2;
#endif

#if ("`$s32L'"==0) && ("`$s32R'"==0)
    id `a'= `b' - 4*i3en + 2 - 2*ctw^2;
#endif

#if ("`$s33L'"==0) && ("`$s33R'"==0)
    id `a'= `b' - 2*i3el - 2 + 2*ctw^2;
#endif

#if ("`$s34L'"==0) && ("`$s34R'"==0)
    id `a'= `b' - 2*i3en + 2 - 2*ctw^2;
#endif

#if ("`$s35L'"==0) && ("`$s35R'"==0)
    id `a'= `b'/stw;
#endif
#if ("`$s36L'"==0) && ("`$s36R'"==0)
    id `a'= `b'*stw;
#endif
*--------------------------------
#if ("`$s37L'"==0)  && ("`$s37R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
#if ("`$s237L'"==0) && ("`$s237R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
#if ("`$s337L'"==0) && ("`$s337R'"==0)
    id `a'= `b'+1;
    argument;
       id `a'= `b'+1;
    endargument;
#endif
*--------------------------------
#if ("`$s38L'"==0)  && ("`$s38R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
#if ("`$s238L'"==0) && ("`$s238R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
#if ("`$s338L'"==0) && ("`$s338R'"==0)
    id `a'= `b'-1;
    argument;
       id `a'= `b'-1;
    endargument;
#endif
*--------------------------------
#if ("`$s39L'"==0)  && ("`$s39R'"==0)
    id `a'= `b' - i3up;
#endif
#if ("`$s239L'"==0) && ("`$s239R'"==0)
    id `a'= `b' - i3ch;
#endif
#if ("`$s339L'"==0) && ("`$s339R'"==0)
    id `a'= `b' - i3tp;
#endif
*--------------------------------
#if ("`$s40L'"==0)  && ("`$s40R'"==0)
    id `a'= `b' + i3up;
#endif
#if ("`$s240L'"==0) && ("`$s240R'"==0)
    id `a'= `b' + i3ch;
#endif
#if ("`$s340L'"==0) && ("`$s340R'"==0)
    id `a'= `b' + i3tp;
#endif
*--------------------------------
#if ("`$s41L'"==0)  && ("`$s41R'"==0)
    id `a'= `b' - i3dn;
#endif
#if ("`$s241L'"==0) && ("`$s241R'"==0)
    id `a'= `b' - i3st;
#endif
#if ("`$s341L'"==0) && ("`$s341R'"==0)
    id `a'= `b' - i3bt;
#endif
*--------------------------------
#if ("`$s42L'"==0)  && ("`$s42R'"==0)
    id `a'= `b' + i3dn;
#endif
#if ("`$s242L'"==0) && ("`$s242R'"==0)
    id `a'= `b' + i3st;
#endif
#if ("`$s342L'"==0) && ("`$s342R'"==0)
    id `a'= `b' + i3bt;
#endif
*--------------------------------
#if ("`$s43L'"==0)  && ("`$s43R'"==0)
    id `a'= `b' - i3el;
#endif
#if ("`$s243L'"==0) && ("`$s243R'"==0)
    id `a'= `b' - i3mo;
#endif
#if ("`$s343L'"==0) && ("`$s343R'"==0)
    id `a'= `b' - i3ta;
#endif
*--------------------------------
#if ("`$s44L'"==0)  && ("`$s44R'"==0)
    id `a'= `b' + i3el;
#endif
#if ("`$s244L'"==0) && ("`$s244R'"==0)
    id `a'= `b' + i3mo;
#endif
#if ("`$s344L'"==0) && ("`$s344R'"==0)
    id `a'= `b' + i3ta;
#endif
*--------------------------------
*45 unused

#if ("`$s46L'"==0) && ("`$s46R'"==0)
    id `a'= 'b' + i3en;
#endif
#endprocedure
*------------
