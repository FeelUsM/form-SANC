#procedure Trace(expressionIn,expressionOut,fourDimension,cIndex)
*****************************************************************
* This procedure evaluates the trace out of the squared matrix  *
* element stored in 'expressionIn'. The summation by the vector *
* boson spin idices are also done;                              *
*                                                               *      
*  expressionIn  -- name of expression to be Traced             *
*  expressionOut -- name of Global expression after Trace       *
*  fourDimension == 1 -- for 4-dimensional Trace                *
*  cIndex -- converts variables into internal notation          *
*****************************************************************

 #call Stop()

 Global TEMP = `expressionIn';            * temporary global for substitutions ...*

.sort :trace-start;
 
 skip;
 nskip TEMP;                 
*==============================	
* -- Step 1: Проверяем, присутствуют ли калибровочно-зависимые термы ...
*============================== 
 #ifdef `debug' 
    #write " ## >>>>>> Checking for a gauge dependence...";
 #endif 
*!!! if (match(xi?sGP) || match(den(xi?sGP,?a)));
*!!!    id xi?sGP$gp = xi;  * $gp is not defined, if 'xi' isn't inside function ('den' here) 
*!!!    print " ## ERROR: there is a gauge parameter \"%$\" in expression: `expressionIn'",$gp; 
*!!!    $error = 1;
*!!!    exit;
*!!! endif;

 #call Stop();
 skip;
 nskip TEMP;
 #call SetFlags()
 
*==============================	
* -- Step 2: используем соотношение полноты для биспиноров
*            и конструируем следы
*==============================            
*
*   -- Step 2-1 : перегруппировываем токи в соотношении TEMP таким образом, что 
*      правильные токи помещаются друг за другом 
*      для использования соотношения полноты и конструирования следов.
*------------------------------
 #ifdef `debug' 
    #write " ## >>>>>> Preparing expression for a trace...";
 #endif 
 #$spIMax = 0;
 if ( match(Vb?suf(?a)) );
     $spI = 1;                 
     repeat;            
* начинаем с самого левого тока Vb(.)*...*U(.)*...*Ub(.)*...*V(.)
       id once        Vb?suf(i?,p?sMom,h?sHel,?a)               =             axgd*Vb(i,p,h,?a);
       repeat id axgd*Vb?suf(i?,p?sMom,h?sHel,?a)*gd?sgd(i?,?b) = gd($spI,?b)*axgd*Vb(i,p,h,?a);
* если начало текущего тока совпало с его концом, создаем след and change $spI for next trace
       if ( match(axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)* V?spf[ppp](i?,p?sMom,h?sHel,?b)) );
           id     axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)* V?spf[ppp](i?,p?sMom,h?sHel,?b) = cR($spI,ppp,p);
           $spI = $spI + 1;                                   
* иначе пытаемся найти другой ток, который следовало бы поставить мосле текущего...
       else;               
           id     axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)*V?spf[ppp1](i?,p1?sMom,h1?sHel,?b) = 
                                                sf(1,ppp,i,p,h)*V(-1,i,p1,h1,?b)*sf(2,ppp1,i,p1,h1);
       endif;
* ...идем его искать
       while ( match(sf(?a)) );
* если присутствует подходящий Vb(..) term...
            if ( match(sf(2,ppp?,i?,p?,h?)*Vb?suf[ppp](j?,p?,h?,?a)) );
* ...then, засовываем term-ы "Vb(..)*..*V(..)" внутрь функции sf(2,..)
                id once sf(2,ppp?,i?,p?,h?)*Vb?suf[ppp](j?,p?,h?,?a) = sf(2,Vb(j,p,h),ppp,i,p,h)*axgd(j);
                repeat id sf(2,?b)*axgd(j?)*gd?sgd(j?,?c)                = sf(2,gd($spI,?c),?b)*axgd(j);
                id once   sf(2,?b)*axgd(j?)*V?spf(j?,p1?sMom,h1?sHel,?c) = sf(2,V(j,p1,h1),?b);
* иначе exit, потому что выражение выглядит некорректным
            else;
               id sf(2,ppp?$x1,?a$x2) = sf(2,ppp,?a);
               $pf = spf[$x1]; 
               $error = 1;  
               print " ## ERROR: can\'t find a pair for  \"%$(%$)\" in expression: `expressionIn'",$pf,$x2;      
               exit;
            endif;
* извлекаем засунутое в sf(2,V(),..,Vb(),?a) после первого тока
            id once   sf(2,U?spf(?a), ?b)*V?spf(-1,?c) = V(-1,?c)*sf(2,?b)*U(-2,?a);
            repeat id sf(2,gd?sgd(?a),?b)*V?spf(-1,?c) = V(-1,?c)*sf(2,?b)*gd(?a);
            id once   sf(2,Vb?suf(?a),?b)*V?spf[ppp](-1,i?,p?,?c) = cR($spI,ppp,p); 
* если конец добавленного тока совпадает с началом первого, создаем след and change $spI for next trace*
            if ( match(sf(1,ppp?,i?,p?,h?)*U?spf[ppp](-2,j?,p?,h?,?a))); 
               id once sf(1,ppp?,i?,p?,h?)*U?spf[ppp](-2,j?,p?,h?,?a) = cR($spI,ppp,p);
               $spI = $spI + 1;             
* иначе продолжаем цикл "while (match(sf(?a)));..."*
            else;
               id once U?spf[ppp](-2,j?,p?,h?,?a) = U(-1,j,p,h,?a)*sf(2,ppp,j,p,h);
            endif;
       endwhile;
     endrepeat;
* сохраняем максимальное значение спиновых индексов. Это потребуется ниже для trace4
     if ($spI > $spIMax) $spIMax = $spI;
 endif;

 #call Stop()
 skip;	
 nskip TEMP;                 
*------------------------------
*   -- Step 2-2 : convert to FORM built in gamma matrices 
*      and use completeness relation for bispinors     
*------------------------------
 id gd(i?,mu?) = g_(i,mu);
 id gd5(i?)    = g5_(i);
 id gd6(i?)    = g6_(i);
 id gd7(i?)    = g7_(i); 
#if `PeskinNaumov'
 id cR(i?,ppp?,p?sMom[mmm]) = g_(i,p) + (-1)^ppp*sMas[mmm]*gI_(i);
#else
 id cR(i?,ppp?,p?sMom[mmm]) = -i_*g_(i,p) + (-1)^ppp*sMas[mmm]*gI_(i);
#endif
 
*==============================
* -- Step 3: use completeness relation for polarization vectors
*============================== 
#if `PeskinNaumov'
**  id pV(23,mu1?,p?sMom,h?sHel,gi1?sAgi)*pVc(23,nu1?,p?sMom,h?sHel,gi2?sAgi) = d_(mu1,nu1); 
  id pV(1,mu1?,p?sMom,h?sHel,?a)*pVc(1,nu1?,p?sMom,h?sHel,?b) = - d_(mu1,nu1);
  id pV(n1?{2,3,-3},mu1?,p?sMom[mmm],h?sHel,?a)*pVc(n1?{2,3,-3},nu1?,p?sMom[mmm],h?sHel,?b) =
                                                     - d_(mu1,nu1) + p(mu1)*p(nu1)/sMas[mmm]/sMas[mmm];
#else
  id pV(23,mu1?,p?sMom,h?sHel,gi1?sAgi)*pVc(23,nu1?,p?sMom,h?sHel,gi2?sAgi) = d_(mu1,nu1); 
  id pV(1,mu1?,p?sMom,h?sHel,?a)*pVc(1,nu1?,p?sMom,h?sHel,?b) = d_(mu1,nu1);
  id pV(n1?{2,3,-3},mu1?,p?sMom[mmm],h?sHel,?a)*pVc(n1?{2,3,-3},nu1?,p?sMom[mmm],h?sHel,?b) =
                                                      d_(mu1,nu1) + p(mu1)*p(nu1)/sMas[mmm]/sMas[mmm];
#endif
						      
* -- Convert masses into the internal notation ...*
 #if ((`cIndex' <= 3) && (`cIndex' >= 0))
     #call Convert(`cIndex')
 #endif

*==============================
* -- Step 4: deal with QCD gauge group generators if any
*==============================
** #call QCDAlgebra();

bracket volum,e,den;
print +s;
.sort :spin-summation-done;
#write "=== SPIN-SUMMED ==="

 skip;
 nskip TEMP; 
*==============================
* -- Step 5: make trace
*==============================
 #ifdef `debug' 
    #write " ## >>>>>> Make trace...";
 #endif 
 if ((match(gd?sgd(?a)) != 0) || (match(Vb?suf(?b)) != 0) || (match(U?spf(?c)) != 0));
    print " ## ERROR in Trace: bad form of expression `expressionIn' :";    
    print " ##       missing bispinor(U,Ub,V,Vb) or not bulitin gamma matrix"; 
    $error = 1;
    exit;
 endif;

 #do i = 1,`$spIMax'-1
    trace4, `i';
    contract;
 #enddo

* go to the 4-dimension ...*
 #if (`fourDimension' == 1)
     id n = 4;
 #endif

* -- Convert masses into the internal notation ...*
 #if ((`cIndex' <= 3) && (`cIndex' >= 0))
     #call Convert(`cIndex')
 #endif
 
 #if (`$isHelicityUsed' == 0)
     #call SetFlags()
 #endif
  
.sort :trace-end-start;
 #call Stop()
 drop TEMP;
 G `expressionOut' = TEMP;
.sort :trace-end-end;

#endprocedure
*------------
