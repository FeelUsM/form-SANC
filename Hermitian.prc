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
#include Misc.h
Cfun	dfun;
sym	ggg;

 skip;
 nskip TEMP;                         * -- skip all the active expresions exept TEMP ...*
*==============================
* -- Step 1: проверяем, присутствуют ли повторяющиеся спиновые индексы
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
* -- Step 2: проверяем на наличие различных спиновых индексов в одном токе (также как в DiracEquation)
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
* -- Step 3: делаем сопряжение векторов поляризации (бозонов)
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
* -- Step 4: сопрягаем биспиноры и матрицы Дирака
*==============================
* Dirac string has a form: Ub?suf(?a)*gd(?c1)*...*gd(?c2)*U?spf(?b)
*
* Notations used:
*   axgd(-1,ppp,?a) --> spf[ppp](?a)	* Set spf:V,U;
*   axgd(-2,uuu,?a) --> suf[uuu](?a)	* Set suf:Vb,Ub;
*   axgd(-3,ggg,?a) --> sgd[ggg](?a)	* Set sgd:gd,gd5,gd6,gd7;
*
 skip;
 nskip TEMP;
 #call GammaRight()

 if ( match(Ub?suf(?a)) );
    repeat;
      id once Ub?suf(?a)*gd?sgd[ggg](?b) = Ub(?a)*axgd(-3,ggg,?b);		* Ub*gd   /once/-> Ub*[gd]
      repeat;
        id axgd(-3,?a)*gd?sgd(?b) = gd(?b)*axgd(-3,?a);				*   [gd_1]*gd_2 -> gd_2*[gd_1]
        id axgd(-3,?a)*U?spf(?b) = U(?b)*axgd(-3,?a);				*   [gd]*U      -> U*[gd]
      endrepeat;
      id Ub?suf[uuu](?a)*U?spf[mmm](?b) = axgd(-2,mmm,?b)*axgd(-1,uuu,?a);	* Ub_1*U_2      -> [Ub_2]*[U_1]
      repeat ;
       #if `PeskinNaumov'
        id axgd(-1,?a)*axgd(-3,1,?b) = + axgd(-3,1,?b)*axgd(-1,?a);		*   [U]*[gd_mu] -> + [gd_mu]*[U]
       #else
        id axgd(-1,?a)*axgd(-3,1,?b) = - axgd(-3,1,?b)*axgd(-1,?a);		*   [U]*[gd_mu] -> - [gd_mu]*[U]
       #endif
        id axgd(-1,?a)*axgd(-3,2,?b) = - axgd(-3,2,?b)*axgd(-1,?a);		*   [U]*[gd5]   -> - [gd5]*[U]
        id axgd(-1,?a)*axgd(-3,3,?b) =   axgd(-3,4,?b)*axgd(-1,?a);		*   [U]*[gd6]   ->   [gd7]*[U]
        id axgd(-1,?a)*axgd(-3,4,?b) =   axgd(-3,3,?b)*axgd(-1,?a);		*   [U]*[gd7]   ->   [gd6]*[U]
      endrepeat; 
    endrepeat;
 endif;

 id axgd(-1,uuu?,?a) = spf[uuu](?a);					* [U]  -> U
 id axgd(-2,uuu?,?a) = suf[uuu](?a);					* [Ub] -> Ub
 id axgd(-3,uuu?,?a) = sgd[uuu](?a);					* [gd] -> gd

*==============================
* -- Step 5: генераторы калибровочной группы QCD
*==============================
** id Tf(gi1?sAgi,cl1?sAcl,cl2?sAcl) = Tf(gi1,cl2,cl1);

*==============================
* -- Step 6: обрабатываем форм-факторы, если есть
*==============================
 #ifdef `formFactorsNumber'
    #do i = 1,`formFactorsNumber'
       #$tmp = sFormFactors[`i'];
       id `$tmp'(?a) = `$tmp'c(?a);
    #enddo
 #endif

*==============================
* -- Step 7: и наконец ...
*==============================
 id i_ = -i_;
** id chi(?a) = chic(?a);
** id den(1,m1?{mw,mz,mh,mtp},Q?) = denc(1,m1,Q);
** id prop(?a) = propc(?a);

.sort
 drop TEMP;
 Global `expressionOut' = TEMP;
.sort

#endprocedure
*------------
