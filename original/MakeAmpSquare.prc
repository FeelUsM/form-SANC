#procedure MakeAmpSquare(expressionIn,expressionOut,multiplyWt)
***************************************************************
* Procedure squares the amplitude "expressionIn", and         *
* stores result in the "expressionOut".                       *
*                                                             *
* "expressionOut" is multiplied by 'multiplyWt'               *
*                                                             *
* NOTE:                                                       *
*   Summation indices (color/lorentz) are treated correctly   *
***************************************************************

 #call Stop()

 skip;
 nskip `expressionIn';

*==============================
* -- Step 1: check whether there are gauge dependent terms
*==============================
*!!! if (match(xi?sGP) || match(den(xi?sGP$gp,?a)));
*!!!    id xi?sGP$gp = xi;   *Note: $gp is not defined, if 'xi' isn't inside function ('den' here)
*!!!    print " ## ERROR: there is a gauge parameter \"%$\" in expression: `expressionIn'",$gp; 
*!!!    $error = 1;
*!!!    exit;
*!!! endif;

*=============================================	
* -- Step 1: make hermitian conjugation. Output is saved in TEMP
*=============================================	
 #call Hermitian(`expressionIn',TEMP);

*=============================================
* -- Step 2: change the Lorentz summation indices
*=============================================
.sort

*// DEBUG - DEBUG - DEBUG - DEBUG
 #$sumIndexMax = 0;
 #$sumIndex = 0;
 #ifndef `$isMakeAmpSquaredUsed'
    #$isMakeAmpSquaredUsed = 0;    * 'true' when this proc. is called first time
 #endif
*// doing this way, we can reset $isMakeAmpSquaredUsed to 0 !
 #if (`$isMakeAmpSquaredUsed' == 0)
     Table relax indexTab(1:1,?a);
      #ifdef `debug'
         Fill indexTab(1) = sf(?a);
      #else
         Fill indexTab(1) = 1;
      #endif   
     #$isMakeAmpSquaredUsed = 1;
    .global
 #endif 
*// DEBUG - DEBUG - DEBUG - DEBUG

 skip;
 nskip TEMP; 

 if ( match(gd(i?,mu?sIndex)) );
     repeat;
       id once gd(i?,mu?sIndex) = axgd(i,mu)*sf(0,mu);       * probe for sum index

       id gd(i?,mu?sIndex[mmm])*sf(0,mu?) = gd(i,sIndexSm[mmm])*sf(1,mu);              * find & change
       id pV?spv(i?,mu?sIndex[mmm],?a)*sf(0,mu?) = pV(i,sIndexSm[mmm],?a)*sf(1,mu);    * find & change
       id pVc?spv(i?,mu?sIndex[mmm],?a)*sf(0,mu?) = pVc(i,sIndexSm[mmm],?a)*sf(1,mu);  * find & change

       id axgd(i?,mu?sIndex[mmm$sumIndex])*sf(1,mu?) = gd(i,sIndexSm[mmm])*indexTab(1,mmm); * if find, change in probe also
       #ifdef `debug'
          if ($sumIndex > $sumIndexMax) $sumIndexMax = $sumIndex;
       #endif
     endrepeat;
     id axgd(i?,mu?sIndex[mmm])*sf(0,mu?) = gd(i,mu);  * not found        
 endif; 
  
 if ( match(pV?spv(i?,mu?sIndex,?a)*pVc?spv(j?,mu?sIndex,?b)) );
     repeat; 
       id once pV?spv(i?,mu?sIndex[mmm$sumIndex],?a)*pVc?spv(j?,mu?sIndex[mmm],?b) = 
                              pV(i,sIndexSm[mmm],?a)*pVc(j,sIndexSm[mmm],?b)*indexTab(1,mmm);
       #ifdef `debug'
          if ($sumIndex > $sumIndexMax) $sumIndexMax = $sumIndex;
       #endif
     endrepeat;
 endif;

*// DEBUG - DEBUG - DEBUG - DEBUG
 #ifdef `debug'
   .sort
    skip;
    nskip TEMP;
    #do i=1,`$sumIndexMax'
       #$index`i' = -1;
       if ( match(sf(`i')) && ($index`i' == -1) );  
           $index`i' = `i';
           $pos = `i';
           $from = sIndex[$pos];
           $to = sIndexSm[$pos];
           print " ## >>>>>> Lorentz Sum Index: %$ --> %$",$from,$to;
       endif;
       id sf(`i') = 1;
    #enddo 
 #endif 
*// DEBUG - DEBUG - DEBUG - DEBUG

*=============================================	
* -- Step 3: change the color summation indices
*=============================================
.sort 

 #$giSMax = 0;
 #$giS = 0;
 #$clSMax = 0;
 #$clS = 0;
 skip;
 nskip TEMP;
 #call QCDAlgebra()
 repeat;
   id once Ta(?a,gi1?sgi[uuu$giS],?b)*Ta(?c,gi1?sgi[uuu],?d) = Ta(?a,scgi[uuu],?b)*Ta(?c,scgi[uuu],?d)*indexTab(1,-1,uuu);
      if ($giS > $giSMax) $giSMax = $giS;
   id once Tf(gi1?sgi[uuu$giS],?a)*Tf(gi1?sgi[uuu],?b) = Tf(scgi[uuu],?a)*Tf(scgi[uuu],?b)*indexTab(1,-1,uuu);
      if ($giS > $giSMax) $giSMax = $giS;
   id once Tf(gi1?sAgi,?a,cl1?scl[uuu$clS],?b)*Tf(gi2?sAgi,?c,cl1?scl[uuu],?d) = Tf(gi1,?a,sccl[uuu],?b)*Tf(gi2,?c,sccl[uuu],?d)*indexTab(1,-2,uuu);
      if ($clS > $clSMax) $clSMax = $clS;
   id Tf(gi1?sgi[uuu$giS],?a)*Ta(?b,gi1?sgi[uuu],?c) = Tf(scgi[uuu],?a)*Ta(?b,scgi[uuu],?c)*indexTab(1,-1,uuu);
      if ($giS > $giSMax) $giSMax = $giS;
 endrepeat; 

*// DEBUG - DEBUG - DEBUG - DEBUG
 #ifdef `debug'
    .sort
     skip;
     nskip TEMP;
     #do i=1,`$giSMax'
        #$indexGI`i' = -1;
        if ( match(sf(-1,`i')) && ($indexGI`i' == -1) );  
            $indexGI`i' = `i';
            print " ## >>>>>> Gluon Color Sum Index: gi`i' --> cgi`i'";
        endif;
        id sf(-1,`i') = 1;
     #enddo 

     #do i=1,`$clSMax'
        #$indexCL`i' = -1;
        if ( match(sf(-2,`i')) && ($indexCL`i' == -1) );  
            $indexCL`i' = `i';
            print " ## >>>>>> Quark Color Sum Index: cl`i' --> ccl`i'";
        endif;
        id sf(-2,`i') = 1;
     #enddo 
 #endif 
*// DEBUG - DEBUG - DEBUG - DEBUG

.sort
 drop TEMP;
*=============================================	
* -- Step 4: construct the square of amplitude `expressionIn'
*=============================================	
 Global `expressionOut' = `multiplyWt'*`expressionIn'*TEMP;
* Global `expressionOut' = TEMP;
 
.sort
 
#endprocedure
*------------
