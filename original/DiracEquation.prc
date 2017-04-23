#procedure DiracEquation(exprName,convertIndex)
***********************************************
*  Uses Dirac equation.                       *
*  exprName     : name of input expression    *
*  convertIndex : converts variables          * 
*                 into internal notation      *
***********************************************

 #call Stop()
 skip;
 nskip `exprName';
*==============================
* -- Step 1: Check whether there are different spinor indices in given current...
*==============================
 id once Vb?suf(ii?sSLI,?a) = Vb(ii,?a)*axgd(ii);
 repeat;
    if (match(axgd(ii?sSLI$x1)*gd?sgd(j?$x2,?b))); 
       id once axgd(ii?sSLI)*gd?sgd(jj?sSLI,?b) = gd(jj,?b)*axgd(ii);
       $xx = theta_(-delta_($x1,$x2),-1);   * = 0 if $x1!=$x2 otherwise = 1 ..*
       if ($xx == 0);
          print " ## ERROR in expression `exprName': %$ != %$ ...",$x1,$x2;
          print " ##        different spin-line index in current..";
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(ii?sSLI$x1)*V?spf(j?$x2,?b)));
       id once axgd(ii?sSLI)*V?spf(jj?sSLI,?b) = V(jj,?b)*axgd(ii);
       $xx = theta_(-delta_($x1,$x2),-1);
       if ($xx == 0);
          print " ## ERROR in expression `exprName': %$ != %$ ...",$x1,$x2;
          print " ##        different spin-line index in current..";
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(ii?sSLI)*Vb?suf(jj?sSLI,?b)));
       id once axgd(ii?sSLI)*Vb?suf(jj?sSLI,?b) = Vb(jj,?b)*axgd(jj);
    endif;     
 endrepeat;
 id axgd(?a) = 1; 

 #call Stop()
 skip;
 nskip `exprName';
 #call SetFlags()
 
 if (match(Vb?suf(ii?sSLI,?a)));
*==============================
* -- Step 2: Introduce an auxilary sf(?a) (CF)function ...
*==============================
    id Vb?suf(ii?sSLI,p?sMom,?a) = Vb(ii,p,?a)*sf(ii,p,?a);

**==============================
* -- Step 3: Check if there are duplicated spinor indecies ...
*==============================
    if (match(sf(ii?sSLI,?a)*sf(ii?sSLI$indx,?b)));
       print " ## ERROR in expression `exprName': %$ - duplicated spinor index ...",$indx;
       print " ##       there are two currents with the same spin-line index...";
       $error = 1;   
       exit;
    endif;

*==============================
* -- Step 4: Invoke GammaRight if there is at least one momentum for Dirac equation ...
*==============================
    if (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));
       #call GammaRight();
       id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];
    endif;

*==============================
* -- Step 5: Use Dirac equation to eliminate corresponding momentum ...
*==============================
    while (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));
      repeat;
*       id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm])*sf(ii?sSLI,p?sMom[mmm],?a) = - sMas[mmm]*sMas[mmm];
        id gd(ii?sSLI,mu?)*gd(ii?sSLI,p?)*sf(ii?sSLI,p?,?a) = (-gd(ii,p)*gd(ii,mu) + 2*p(mu))*sf(ii,p,?a);
        id Vb?suf[uuu](ii?sSLI,p?sMom[mmm],?a)*gd(ii?sSLI,p?sMom[mmm]) = (-1)^uuu*i_*sMas[mmm]*Vb(ii,p,?a);
        id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];
      endrepeat;
    endwhile;

*==============================
* -- Step 6: Kick out the auxilary sf(?a) function ...
*==============================
    id Vb?suf(ii?sSLI,p?sMom,?a)*sf(ii?sSLI,p?sMom,?a) = Vb(ii,p,?a);

 endif;

* -- Convert appeared masses to the internal notation ...*
 #if ((`convertIndex' <= 3) && (`convertIndex' >= 0))
     #call Convert(`convertIndex')
 #endif
 #call SetFlags() 

.sort

*==============================
* -- Step 7: Do the same procedure for the particle field(spf) ...
*==============================
 #call Stop()
 skip;
 nskip `exprName';
 if (match(U?spf(ii?sSLI,?a)));

    id U?spf(ii?sSLI,p?sMom,?a) = U(ii,p,?a)*sf(ii,p,?a);

    if (match(sf(ii?sSLI,?a)*sf(ii?sSLI$indx,?b)));
       print " ## ERROR in expression `exprName': %$ - duplicated spinor index ...",$indx;
       print " ##       there are two currents with the same spin-line index...";
       $error = 1;
       exit;
    endif;

    if (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));
       #call GammaLeft()
       id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];
    endif;

    while (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));
       repeat;
*        id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm])*sf(ii?sSLI,p?sMom[mmm],?a) = - sMas[mmm]*sMas[mmm];
         id gd(ii?sSLI,p?)*gd(ii?sSLI,mu?)*sf(ii?sSLI,p?,?a) = (-gd(ii,mu)*gd(ii,p) + 2*p(mu))*sf(ii,p,?a);
         id gd(ii?sSLI,p?sMom[mmm])*U?spf[ppp](ii?sSLI,p?sMom[mmm],?a) = (-1)^ppp*i_*sMas[mmm]*U(ii,p,?a);
         id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];
       endrepeat;
    endwhile;

    id U?spf(ii?sSLI,p?sMom,?a)*sf(ii?sSLI,p?sMom,?a) = U(ii,p,?a);

 endif;

* -- convert appeared masses to the internal notation ...*
 #if ((`convertIndex' <= 3) && (`convertIndex' >= 0))
     #call Convert(`convertIndex')
 #endif
 #call GammaRight();
 #call SetFlags() 
 #call Stop()

#endprocedure
*------------
