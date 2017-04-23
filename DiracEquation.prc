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
* -- Step 1: проверка на наличие разных спиновых индексов в одном токе...
*==============================
 id once Vb?suf(ii?sSLI,?a) = Vb(ii,?a)*axgd(ii);			* Vb(ii) /once/-> Vb(ii)*axgd(ii)
 repeat;
    if (match(axgd(ii?sSLI$x1)*gd?sgd(j?$x2,?b)));			* if_match axgd(ii)*gd(j)
       id once axgd(ii?sSLI)*gd?sgd(jj?sSLI,?b) = gd(jj,?b)*axgd(ii);	* axgd(ii)*gd(j) /once/-> gd(j)*axgd(ii)
       $xx = theta_(-delta_($x1,$x2),-1);   				* = 0 if $x1!=$x2 otherwise = 1 ..*
       if ($xx == 0);
          print " ## ERROR in expression `exprName': %$ != %$ ...",$x1,$x2;
          print " ##        different spin-line index in current..";
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(ii?sSLI$x1)*V?spf(j?$x2,?b)));			* if_match axgd(ii)*V(j)
       id once axgd(ii?sSLI)*V?spf(jj?sSLI,?b) = V(jj,?b)*axgd(ii);	* axgd(ii)*V(jj) /once/-> V(jj)*axgd(ii)
       $xx = theta_(-delta_($x1,$x2),-1);
       if ($xx == 0);
          print " ## ERROR in expression `exprName': %$ != %$ ...",$x1,$x2;
          print " ##        different spin-line index in current..";
          $error = 1;
          exit;
       endif;
    endif;
    if (match(axgd(ii?sSLI)*Vb?suf(jj?sSLI,?b)));			* if_match axgd(ii)*Vb(j)
       id once axgd(ii?sSLI)*Vb?suf(jj?sSLI,?b) = Vb(jj,?b)*axgd(jj);	* axgd(ii)*Vb(j) /once/->Vb(j)*axgd(j)
    endif;     
 endrepeat;
 id axgd(?a) = 1; 							* axgd() -> 1

 #call Stop()
 skip;
 nskip `exprName';
 #call SetFlags()
 
 if (match(Vb?suf(ii?sSLI,?a)));					* по крайней мере одна свертка гамма-матриц/спиноров
*==============================
* -- Step 2: вводим дополнительную коммутирующую функцию sf(?a) (CF)function ...
*==============================
    id Vb?suf(ii?sSLI,p?sMom,?a) = Vb(ii,p,?a)*sf(ii,p,?a);		* Vb(ii) -> Vb(ii,p)*sf(ii,p)

**==============================
* -- Step 3: ищем повторяющиеся спиновые индексы ...
*==============================
    if (match(sf(ii?sSLI,?a)*sf(ii?sSLI$indx,?b)));			* if_match sf(ii)*sf(ii)
       print " ## ERROR in expression `exprName': %$ - duplicated spinor index ...",$indx;
       print " ##       there are two currents with the same spin-line index...";
       $error = 1;   
       exit;
    endif;

*==============================
* -- Step 4: Вызываем GammaRight если присутствует по крайней мере один момент с крышкой ...
*==============================
    if (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));				* if_match sf(ii,p)*gd(ii,p)
       #call GammaRight();								* GammaRight();
       id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];	* \hat p*\hat p -> -m^2
    endif;

*==============================
* -- Step 5: Используем уравнение Дирака, чтобы выкинуть соответствующие моменты ...
*==============================
    while (match(sf(ii?sSLI,p?sMom,?a)*gd(ii?sSLI,p?sMom)));
      repeat;
*       id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm])*sf(ii?sSLI,p?sMom[mmm],?a) = - sMas[mmm]*sMas[mmm];
        id gd(ii?sSLI,mu?)*gd(ii?sSLI,p?)*sf(ii?sSLI,p?,?a) = (-gd(ii,p)*gd(ii,mu) + 2*p(mu))*sf(ii,p,?a);
        id Vb?suf[uuu](ii?sSLI,p?sMom[mmm],?a)*gd(ii?sSLI,p?sMom[mmm]) = (-1)^uuu*i_*sMas[mmm]*Vb(ii,p,?a);
        id gd(ii?sSLI,p?sMom[mmm])*gd(ii?sSLI,p?sMom[mmm]) = - sMas[mmm]*sMas[mmm];
      endrepeat;
    endwhile;
*    while_match sf(ii,p)*gd(ii,p);
*      repeat;
**       gd(p)*gd(p)*sf(p)  -> - m^2;
*        gd(mu)*gd(p)*sf(p) -> (-gd(p)*gd(mu) + 2*p(mu))*sf(p);
*        Vb(p)*gd(p)        -> (+/- = Ub/Vb)*i_*m*Vb(p);
*        gd(p)*gd(p)        -> - m^2;
*      endrepeat;
*    endwhile;

*==============================
* -- Step 6: Выкидываеи дополнительную функцию sf(?a) ...
*==============================
    id Vb?suf(ii?sSLI,p?sMom,?a)*sf(ii?sSLI,p?sMom,?a) = Vb(ii,p,?a);		* Vb(ii,p)*sf(ii,p) -> Vb(ii)

 endif;

* -- Конвертируем полученные массы во внутреннюю нотацию ...*
 #if ((`convertIndex' <= 3) && (`convertIndex' >= 0))
     #call Convert(`convertIndex')
 #endif
 #call SetFlags() 

.sort

*==============================
* -- Step 7: Делаем ту же процедуру для "правых" полей(spf) ...
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
