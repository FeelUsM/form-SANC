#procedure QCDAlgebra()
********************************************************
* T^{a} (a=1,...,8) --  SU(3) color-group generators;  *
*                                                      *
* Properties:                                          *
* ===========                                          *
*  T^{+} = T, Tr( T^{a} ) = 0                          *
*                                                      *
*  [T^{a},T^{b}] = i_f^{abc}T^{c},                     *
*  f^{abc} -- structural constants of SU(3)            *
*                                                      *
*  Tr(T^{a} T^{b}) = d_(a,b)/2                         *
*  sum_{a}( T^{a} T^{a} ) = Cf                         *
*                                                      *   
*  sum_{a,b}( f^{abc} f^{abd} ) = Ca d_(c,d)           *
*                                                      *
*  sum_{a}( T^{a} T^{b} T^{a} ) = (Cf - Ca/2) T^{b}    *
*  sum_{b,c}( f^{abc}T^{b}T^{c} ) = i_/2 Ca T^{a}      *
*                                                      *
*  Ca = Nc, Cf = (Nc^2-1)/2/Nc;                        *
*                                                      *
* Notations:                                           *
* ----------                                           *
*  Tf == T            ->  fundamental representation   *
*  Ta == -i_*f^{abc}  ->  adjoint representation       *
*  Nc = 3             ->  QCD color group dimension    *
********************************************************
 repeat;

*// Tr ( T^{gi1} ) = 0
   id Tf(gi1?sAgi,cl1?sAcl,cl1?sAcl) = 0;

*// t^{gi1}*t^{gi2}*t^{gi1} = (-1/2*Ca + Cf)*t^{gi1}
   id Tf(gi1?sAgi,cl1?sAcl,cl2?sAcl)*Tf(gi2?sAgi,cl2?sAcl,cl3?sAcl)*
                                     Tf(gi1?sAgi,cl3?sAcl,cl4?sAcl) = (Cf - 1/2*Ca)*Tf(gi2,cl1,cl4);

*// f^{abc}*t^{b}*t^{c} = i_/2*Ca*t^{a}
   repeat; 
     id Ta(gi1?sAgi,?a)*Tf(gi1?sAgi,?b) = Ta(?a,gi1)*Tf(gi1,?b) ;
     id Ta(gi1?sAgi,gi2?sAgi,gi3?sAgi)*Tf(gi2?sAgi,cl3?sAcl,cl5?sAcl)*
                                       Tf(gi3?sAgi,cl5?sAcl,cl2?sAcl) = + i_/2*Ca*Tf(gi1,cl3,cl2);
     id Ta(gi1?sAgi,gi3?sAgi,gi2?sAgi)*Tf(gi2?sAgi,cl3?sAcl,cl5?sAcl)*
                                       Tf(gi3?sAgi,cl5?sAcl,cl2?sAcl) = - i_/2*Ca*Tf(gi1,cl3,cl2);
   endrepeat;

*// Tr (t^{gi1}*t^{gi2}) = d_(gi1,gi2)/2
   id Tf(gi1?sAgi,cl2?sAcl,cl3?sAcl)*Tf(gi2?sAgi,cl3?sAcl,cl2?sAcl) = d(gi1,gi2)/2;

*// sum_{gi1} t^{gi1}*t^{gi1) = Cf
   id Tf(gi1?sAgi,cl1?sAcl,cl2?sAcl)*Tf(gi1?sAgi,cl2?sAcl,cl3?sAcl) = d(cl1,cl3)*Cf;

*// f^{abc}f^{abd} = Ca d_(c,d)
   repeat;
     id Ta(gi1?sAgi,?a)*Ta(?b,gi1?sAgi,?c) = Ta(?a,gi1)*Ta(?b,gi1,?c);
     id Ta(gi1?sAgi,gi2?sAgi,gi3?sAgi)*Ta(gi4?sAgi,gi3?sAgi,gi2?sAgi) = -Ta(gi1,gi2,gi3)*Ta(gi4,gi2,gi3);
     id Ta(gi1?sAgi,gi2?sAgi,gi3?sAgi)*Ta(gi4?sAgi,gi2?sAgi,gi3?sAgi) = Ca*d(gi1,gi4);
   endrepeat;

*// sum indices by Cronecker-delta function 
   id Tf(gi1?sAgi,?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = Tf(gi1,?a,?c,?d,?b);
   id Tf(gi1?sAgi,?a)*d(?c,gi1?sAgi,?d) = Tf(?c,?d,?a);
   id d(?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = d(?a,?c,?d,?b);

   id Ta(?a,gi1?sAgi,?b)*d(?c,gi1?sAgi,?d) = Ta(?a,?c,?d,?b);
   id d(?a,gi1?sAgi,?b)*d(?c,gi1?sAgi,?d) = d(?a,?c,?d,?b);

 endrepeat;
  
*// sum by gluon/quark color 
 id d(gi1?sAgi,gi1?sAgi) = 2*Nc*Cf;  
 id d(cl1?sAcl,cl1?sAcl) = Nc;

#endprocedure
*------------
