#procedure Hermitian(expressionIn,expressionOut) 
*************************************************
* This procedure makes hermitian conjugation    *
* of "expressionIn".                            *
*                                               *
* Output is stored in Global "expressionOut"    *
*************************************************

 #call Stop()

 Global TEMP = `expressionIn';       * -- temporary global for substitutions ...*

.sort

 skip;
 nskip TEMP;                         * -- skip all the active expresions exept TEMP ...*
*==============================
* -- Step 1: check whether there are duplicated spinor indices
*==============================
 id Vb?suf(i?,p?sMom,h?sHel,?a) = Vb(i,p,h,?a)*sf(i,p,?a);
 if (match(sf(i?,?a)*sf(i?$indx,?b)));
    print " ## ERROR in expression `expressionIn':";
    print " ##       \"%$\" --> duplicated spin-line index ...",$indx;
    $error = 1;
 else;
    repeat id sf(?a) = 1;
 endif;

 #call Stop()
 skip;
 nskip TEMP;                 
*==============================
* -- Step 2: check whether there are different spinor indices in given current
*==============================
 id once Vb?suf(i?,?a) = Vb(i,?a)*axgd(i);
 repeat; 
    if (match(axgd(i?$x1)*gd?sgd(ii?$x2,?b))); 
       id once axgd(i?)*gd?sgd(ii?,?b) = gd(ii,?b)*axgd(i);
       $xx = theta_(-delta_($x1,$x2),-1);   * = 0 if $x1!=$x2 otherwise = 1 ..*
       if ($xx == 0);
          print " ## ERROR in expression `expressionIn':";
          print "        different spin-line indices in current: %$ != %$ ...",$x1,$x2;
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(i?$x1)*V?spf(ii?$x2,?b))); 
       id once axgd(i?)*V?spf(ii?,?b) = V(ii,?b)*axgd(i);
       $xx = theta_(-delta_($x1,$x2),-1);
       if ($xx == 0);
          print " ## ERROR in expression `expressionIn':";
          print " ##       different spin-line indices in current: %$ != %$ ...",$x1,$x2;
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(i?)*Vb?suf(ii?,?b))); 
       id once axgd(i?)*Vb?suf(ii?,?b) = Vb(ii,?b)*axgd(ii);
    endif;     
 endrepeat;
 id axgd(?a) = 1;

 #call Stop()
*==============================
* -- Step 3: make conjugation of polarization vectors
*==============================
 skip;
 nskip TEMP;
 id pV?spv[mmm](?a) = dfun(mmm,?a);
.sort
 skip;
 nskip TEMP;
 id dfun(1,?a) = spv[2](?a);
 id dfun(2,?a) = spv[1](?a);

.sort
 
*==============================
* -- Step 4: conjugate bispinors and Dirac gamma matrices
*==============================
* Dirac string has a form: Ub?suf(?a)*gd(?c1)*...*gd(?c2)*U?spf(?b)
*
* Notations used:
*   axgd(-1,ppp,?a) --> spf[ppp](?a)
*   axgd(-2,uuu,?a) --> suf[uuu](?a)
*   axgd(-3,ggg,?a) --> sgd[ggg](?a)
*
 skip;
 nskip TEMP;
 #call GammaRight()

 if ( match(Ub?suf(?a)) );
    repeat;
      id once Ub?suf(?a)*gd?sgd[ggg](?b) = Ub(?a)*axgd(-3,ggg,?b);
      repeat;
        id axgd(-3,?a)*gd?sgd(?b) = gd(?b)*axgd(-3,?a);
        id axgd(-3,?a)*U?spf(?b) = U(?b)*axgd(-3,?a);
      endrepeat;
      id Ub?suf[uuu](?a)*U?spf[mmm](?b) = axgd(-2,mmm,?b)*axgd(-1,uuu,?a);
      repeat ;
        id axgd(-1,?a)*axgd(-3,1,?b) = - axgd(-3,1,?b)*axgd(-1,?a);
        id axgd(-1,?a)*axgd(-3,2,?b) = - axgd(-3,2,?b)*axgd(-1,?a);
        id axgd(-1,?a)*axgd(-3,3,?b) =   axgd(-3,4,?b)*axgd(-1,?a);
        id axgd(-1,?a)*axgd(-3,4,?b) =   axgd(-3,3,?b)*axgd(-1,?a);
      endrepeat; 
    endrepeat;
 endif;

 id axgd(-1,uuu?,?a) = spf[uuu](?a);
 id axgd(-2,uuu?,?a) = suf[uuu](?a);
 id axgd(-3,uuu?,?a) = sgd[uuu](?a);

*==============================
* -- Step 5: the QCD gauge group generators
*==============================
 id Tf(gi1?sAgi,cl1?sAcl,cl2?sAcl) = Tf(gi1,cl2,cl1);

*==============================
* -- Step 6: deal with the Form Factors, if any 
*==============================
 #ifdef `formFactorsNumber'
    #do i = 1,`formFactorsNumber'
       #$tmp = sFormFactors[`i'];
       id `$tmp'(?a) = `$tmp'c(?a);
    #enddo
 #endif

*==============================
* -- Step 7: and finally ...
*==============================
 id i_ = -i_;
 id chi(?a) = chic(?a);
 id den(1,m1?{mw,mz,mh,mtp},Q?) = denc(1,m1,Q);
 id prop(?a) = propc(?a);

.sort
 drop TEMP;
 Global `expressionOut' = TEMP;
.sort

#endprocedure
*------------
