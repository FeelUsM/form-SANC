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

.sort
 
 skip;
 nskip TEMP;                 
*==============================	
* -- Step 1: Check whether there are gauge dependent terms ...
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
* -- Step 2: use the completeness relation for bispinors 
*            and construct the traces
*==============================            
*
*   -- Step 2-1 : regroup the currents in the TEMP expression in such 
*      a way that the proper currents are placed side by side in order 
*      to use the completeness relation and construct traces.
*------------------------------
 #ifdef `debug' 
    #write " ## >>>>>> Preparing expression for a trace...";
 #endif 
 #$spIMax = 0;
 if ( match(Vb?suf(?a)) );
     $spI = 1;                 
     repeat;            
* begin from the first left current Vb(.)*...*U(.)*...*Ub(.)*...*V(.)
       id once Vb?suf(i?,p?sMom,h?sHel,?a) = axgd*Vb(i,p,h,?a);
       repeat id axgd*Vb?suf(i?,p?sMom,h?sHel,?a)*gd?sgd(i?,?b) = gd($spI,?b)*axgd*Vb(i,p,h,?a);
* if this completes the whole trace do it and change $spI for next trace
       if ( match(axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)*V?spf[ppp](i?,p?sMom,h?sHel,?b)) );
           id axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)*V?spf[ppp](i?,p?sMom,h?sHel,?b) = cR($spI,ppp,p);
           $spI = $spI + 1;                                   
* esle try to find other currnet, which completes the trace...
       else;               
           id axgd*Vb?suf[ppp](i?,p?sMom,h?sHel,?a)*V?spf[ppp1](i?,p1?sMom,h1?sHel,?b) = 
                                                 sf(1,ppp,i,p,h)*V(-1,i,p1,h1,?b)*sf(2,ppp1,i,p1,h1);
       endif;
* ...so, try to find
       while ( match(sf(?a)) );
* if there is a proper Vb(..) term...
            if ( match(sf(2,ppp?,i?,p?,h?)*Vb?suf[ppp](j?,p?,h?,?a)) );
* ...then, collect terms "Vb(..)*..*V(..)" inside the function sf(2,..)
                id once sf(2,ppp?,i?,p?,h?)*Vb?suf[ppp](j?,p?,h?,?a) = sf(2,Vb(j,p,h),ppp,i,p,h)*axgd(j);
                repeat id sf(2,?b)*axgd(j?)*gd?sgd(j?,?c) = sf(2,gd($spI,?c),?b)*axgd(j);
                id once sf(2,?b)*axgd(j?)*V?spf(j?,p1?sMom,h1?sHel,?c) = sf(2,V(j,p1,h1),?b);
* esle exit, because the expression seems to be incorrect
            else;
               id sf(2,ppp?$x1,?a$x2) = sf(2,ppp,?a);
               $pf = spf[$x1]; 
               $error = 1;  
               print " ## ERROR: can\'t find a pair for  \"%$(%$)\" in expression: `expressionIn'",$pf,$x2;      
               exit;
            endif;
* move the collected sf(2,V(),..,Vb(),?a) to the first current and extract*
            id once sf(2,U?spf(?a),?b)*V?spf(-1,?c) = V(-1,?c)*sf(2,?b)*U(-2,?a);
            repeat id sf(2,gd?sgd(?a),?b)*V?spf(-1,?c) = V(-1,?c)*sf(2,?b)*gd(?a);
            id once sf(2,Vb?suf(?a),?b)*V?spf[ppp](-1,i?,p?,?c) = cR($spI,ppp,p); 
* if this completes the hole trace do it and change $spI for next trace*
            if ( match(sf(1,ppp?,i?,p?,h?)*U?spf[ppp](-2,j?,p?,h?,?a))); 
               id once sf(1,ppp?,i?,p?,h?)*U?spf[ppp](-2,j?,p?,h?,?a) = cR($spI,ppp,p);
               $spI = $spI + 1;             
* else proceed the loop "while (match(sf(?a)));..."*
            else;
               id once U?spf[ppp](-2,j?,p?,h?,?a) = U(-1,j,p,h,?a)*sf(2,ppp,j,p,h);
            endif;
       endwhile;
     endrepeat;
* save the maximal value of the spin line index. Need this below in trace4
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
 id cR(i?,ppp?,p?sMom[mmm]) = -i_*g_(i,p) + (-1)^ppp*sMas[mmm]*gI_(i);
 
*==============================
* -- Step 3: use completeness relation for polarization vectors
*============================== 
  id pV(23,mu1?,p?sMom,h?sHel,gi1?sAgi)*pVc(23,nu1?,p?sMom,h?sHel,gi2?sAgi) = d_(mu1,nu1); 
  id pV(1,mu1?,p?sMom,h?sHel,?a)*pVc(1,nu1?,p?sMom,h?sHel,?b) = d_(mu1,nu1);
  id pV(s?{2,3,-3},mu1?,p?sMom[mmm],h?sHel,?a)*pVc(s?{2,3,-3},nu1?,p?sMom[mmm],h?sHel,?b) =
                                                      d_(mu1,nu1) + p(mu1)*p(nu1)/sMas[mmm]/sMas[mmm];
 
* -- Convert masses into the internal notation ...*
 #if ((`cIndex' <= 3) && (`cIndex' >= 0))
     #call Convert(`cIndex')
 #endif

*==============================
* -- Step 4: deal with QCD gauge group generators if any
*==============================
 #call QCDAlgebra();

.sort

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
  
.sort
 #call Stop()
 drop TEMP;
 G `expressionOut' = TEMP;
.sort

#endprocedure
*------------
