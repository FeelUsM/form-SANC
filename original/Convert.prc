#procedure Convert(index)
*********************/
* This procedure converts variables to  *
* the internal notation.                *
*                                       *
*   index = 0 : convert ALL variables   *
*   index = 1 : convert ONLY masses     *
*   index = 2 : convert Energy and Beta *
*   index = 3 : convert th and ph       *
*                                       *
*********************/
 #if (`$num' != 0)
     #$nLoop = `$num';
 #else
     #$nLoop = `nMomenta';
 #endif
 #do i=1,`$nLoop'

    #ifdef `Im`i'';
       #if ((`index' == 1) || (`index' == 0))
           argument;
              id m`i' = pm(`I`i'');                          * -- convert "m"asses
           endargument;
           id m`i'^n?!{,0} = pm(`I`i'')^n;                         
       #endif
       #$particleName = pn(`I`i'');
       #if ((`index' == 2) || (`index' == 0))
           id beta`i' = beta`$particleName';                 * -- convert "beta"s
           id E`i' = E`$particleName';                       * -- convert "E"nergies
       #endif
       #if ((`index' == 3) || (`index' == 0))
            argument;
              id th`i' = th`$particleName';                  * -- convert "th" angles
              id ph`i' = ph`$particleName';                  * -- convert "ph" angles
            endargument;
       #endif
    #endif
 
    #ifdef `Iin`i'';
       #if ((`index' == 1) || (`index' == 0))
           argument;
              id m`i' = pm(`Iin`i'');                        * -- convert "m"asses
           endargument;
           id m`i'^n?!{,0} = pm(`Iin`i'')^n;                         
       #endif
       #$particleName = pn(`Iin`i'');
       #if ((`index' == 2) || (`index' == 0))
           id beta`i' = beta`$particleName';                 * -- convert "beta"s
           id E`i'    = E`$particleName';                    * -- convert "E"nergies
       #endif
       #if ((`index' == 3) || (`index' == 0))
           argument;
              id th`i' = th`$particleName';                  * -- convert "th" angles
              id ph`i' = ph`$particleName';                  * -- convert "ph" angles
           endargument;
       #endif
    #endif

    #ifdef `Ifn`i'';
       #if ((`index' == 1) || (`index' == 0))
           argument;
              id m`i' = pm(`Ifn`i'');                        * -- convert "m"asses
           endargument;
           id m`i'^n?!{,0} = pm(`Ifn`i'')^n;                         
       #endif
       #$particleName = pn(`Ifn`i'');
       #if ((`index' == 2) || (`index' == 0))
           id beta`i' = beta`$particleName';                 * -- convert "beta"s
           id E`i'    = E`$particleName';                    * -- convert "E"nergies
       #endif
       #if ((`index' == 3) || (`index' == 0))
           argument;
              id th`i' = th`$particleName';                  * -- convert "th" angles
              id ph`i' = ph`$particleName';                  * -- convert "ph" angles
           endargument;
       #endif
    #endif

 #enddo
 id mgm = 0;
 id mgn = 0;
 argument;
    id mgm = 0;
    id mgn = 0;
 endargument;
 id sr2^-2 = 1/2;
 id sr2^2  = 2;
#endprocedure
*------------
