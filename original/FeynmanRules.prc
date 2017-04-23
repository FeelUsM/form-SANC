#procedure FeynmanRules(I)
*----------------------------
* I=0: QED,QCD
* I=1: EW
*
*
* ##########################################################################################
* #                                                                                        #
* #                            FEYNMAN RULES FOR QED,QCD and EW                            #
* #                                                                                        #
* ##########################################################################################
* 
*
* ###################################################
* #          GLUON-GLUON/GHOST VERTICES             #
* ###################################################
*
* ===== Tree Point Vertices =====
* ===============================
*
*             al,gi2,q
*               /
*              /
*            \/_
*   ___\_____/
*      /     \
* mu,gi1,p    \_
*             /\
*               \
*             be,gi3,k
*
id vert(+23,+23,-23,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*(+d_(mu,al)*(p(be)-q(be))
                                                                                         +d_(al,be)*(q(mu)-k(mu))
                                                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(+23,-23,+23,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*(+d_(mu,al)*(p(be)-q(be))
                                                                                         +d_(al,be)*(q(mu)-k(mu))
                                                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(-23,+23,-23,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*(+d_(mu,al)*(p(be)-q(be))
                                                                                         +d_(al,be)*(q(mu)-k(mu))
                                                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(-23,-23,+23,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*(+d_(mu,al)*(p(be)-q(be))
                                                                                         +d_(al,be)*(q(mu)-k(mu))
                                                                                         +d_(mu,be)*(k(al)-p(al)));
*
*
* ===== Four Point Vertices =====
* ===============================
*
*   l      k
* la,gi2  de,gi3
*    \     /
*     \   /
*      \ /
*       X
*      / \
*     /   \
*    /     \
* mu,gi1  nu,gi4
*   i      j
*        l i   k   j   
id vert(23,23,-23,-23,la?,mu?,de?,nu?,gi1?sgi,gi2?sgi,gi3?sgi,gi4?sgi) = -i_*gs^2*(
				   +Ta(gi2,gi3,gi5)*Ta(gi1,gi4,gi5)*(+d_(la,mu)*d_(nu,de)
                                          		             -d_(nu,la)*d_(mu,de))
                                   +Ta(gi2,gi1,gi5)*Ta(gi3,gi4,gi5)*(+d_(la,de)*d_(mu,nu)
                                                                     -d_(nu,la)*d_(mu,de))
				   +Ta(gi2,gi4,gi5)*Ta(gi3,gi1,gi5)*(+d_(la,de)*d_(mu,nu)
                                                                     -d_(mu,la)*d_(de,nu))
                                                      );

id vert(23,23,+23,+23,la?,mu?,de?,nu?,gi2?sgi,gi1?sgi,gi3?sgi,gi4?sgi) = -i_*gs^2*(
				   +Ta(gi2,gi3,gi5)*Ta(gi1,gi4,gi5)*(+d_(la,mu)*d_(nu,de)
                                                                     -d_(nu,la)*d_(mu,de))
                                   +Ta(gi2,gi1,gi5)*Ta(gi3,gi4,gi5)*(+d_(la,de)*d_(mu,nu)
                                                                     -d_(nu,la)*d_(mu,de))
                                   +Ta(gi2,gi4,gi5)*Ta(gi3,gi1,gi5)*(+d_(la,de)*d_(mu,nu)
                                                                     -d_(mu,la)*d_(de,nu))
                                                      );
*
*
* ===== Gluon Ghost Vertices =====
* ================================
*
id vert(+23,+24,-24,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi3,gi2)*k(mu);
id vert(+23,-24,+24,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*q(mu);
id vert(-23,+24,-24,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi3,gi2)*k(mu);
id vert(-23,-24,+24,mu?,al?,be?,gi1?sgi,gi2?sgi,gi3?sgi,p?,q?,k?) = -gs*Ta(gi1,gi2,gi3)*q(mu);
*
*
*
* ###################################################
* #             QUARK-GLUON VERTICES                #
* ###################################################
*
*              cl2
*              /
*             /
*           \/_
*   ________/
* mu,gi1    \
*            \
*            _\/
*              \
*             cl3
*
id vert(n?{23,-23},+13,-13,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-13,+13,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},+14,-14,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-14,+14,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},+17,-17,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-17,+17,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},+18,-18,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-18,+18,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-21,+21,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},+21,-21,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},-22,+22,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
id vert(n?{23,-23},+22,-22,mu?,gi1?sgi,cl2?scl,cl3?scl,ii?) = i_/2*gs*gd(ii,mu)*(gd6(ii)+gd7(ii))*Tf(gi1,cl2,cl3);
*
*
*
* ###########################################################################################
* #                                                                                         #
* #                                FEYNMAN RULES FOR SM                                     #
* #                                                                                         #
* ###########################################################################################
*
*             -------------boson_boson_vertices-----------
*
id vert(n?{1,-1},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{1,-1},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
*
id vert(n?{1,-1},-3,+6,mu?,al?,be?,p?,q?,k?) = -i_*g*stw*pm(3)*d_(mu,al);
id vert(n?{1,-1},+6,-3,mu?,al?,be?,p?,q?,k?) = -i_*g*stw*pm(3)*d_(mu,be);
id vert(n?{1,-1},+3,-6,mu?,al?,be?,p?,q?,k?) = +i_*g*stw*pm(3)*d_(mu,al);
id vert(n?{1,-1},-6,+3,mu?,al?,be?,p?,q?,k?) = +i_*g*stw*pm(3)*d_(mu,be);
*
id vert(n?{1,-1},-6,+6,mu?,al?,be?,p?,q?,k?) =    +g*stw*(k(mu)-q(mu));
id vert(n?{1,-1},+6,-6,mu?,al?,be?,p?,q?,k?) =    +g*stw*(q(mu)-k(mu));
id vert(n?{1,-1},-7,+7,mu?,al?,be?,p?,q?,k?) =    -g*stw*q(mu)/rxi(3);
id vert(n?{1,-1},+7,-7,mu?,al?,be?,p?,q?,k?) =    -g*stw*k(mu)/rxi(3);
id vert(n?{1,-1},-8,+8,mu?,al?,be?,p?,q?,k?) =    +g*stw*q(mu)/rxi(3);
id vert(n?{1,-1},+8,-8,mu?,al?,be?,p?,q?,k?) =    +g*stw*k(mu)/rxi(3);
*
*-
*
id vert(n?{2,-2},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{2,-2},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n1?{2,-2},n2?{5,-5},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(q(mu)-k(mu));
id vert(n1?{2,-2},n2?{4,-4},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(k(mu)-q(mu));
id vert(n1?{2,-2},n2?{2,-2},n3?{4,-4},mu?,al?,be?,p?,q?,k?) =    -g*pm(3)/ctw^2*d_(mu,al);
id vert(n1?{2,-2},n2?{4,-4},n3?{2,-2},mu?,al?,be?,p?,q?,k?) =    -g*pm(3)/ctw^2*d_(mu,be);
*
id vert(n?{2,-2},-3,+6,mu?,al?,be?,p?,q?,k?) = +i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(n?{2,-2},+6,-3,mu?,al?,be?,p?,q?,k?) = +i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(n?{2,-2},-6,+3,mu?,al?,be?,p?,q?,k?) = -i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(n?{2,-2},+3,-6,mu?,al?,be?,p?,q?,k?) = -i_*g*stw^2/ctw*pm(3)*d_(mu,al);
*
id vert(n?{2,-2},-6,+6,mu?,al?,be?,p?,q?,k?) =    +g*(2*ctw-1/ctw)/2*(k(mu)-q(mu));
id vert(n?{2,-2},+6,-6,mu?,al?,be?,p?,q?,k?) =    +g*(2*ctw-1/ctw)/2*(q(mu)-k(mu));
id vert(n?{2,-2},-7,+7,mu?,al?,be?,p?,q?,k?) =    -g*ctw*q(mu)/rxi(3);
id vert(n?{2,-2},+7,-7,mu?,al?,be?,p?,q?,k?) =    -g*ctw*k(mu)/rxi(3);
id vert(n?{2,-2},-8,+8,mu?,al?,be?,p?,q?,k?) =    +g*ctw*q(mu)/rxi(3);
id vert(n?{2,-2},+8,-8,mu?,al?,be?,p?,q?,k?) =    +g*ctw*k(mu)/rxi(3);
*
*-
*
* ((1,2)33)--(pqk) ---> (kpq)--(33(1,2))
id vert(+3,-3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,+3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
* ((1,2)33)--(pqk) ---> (qkp)--(3(1,2)3)
id vert(+3,n?{1,-1},-3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))    
                                                          +d_(mu,be)*(k(al)-p(al))); 
id vert(-3,n?{1,-1},+3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))      
                                                          +d_(mu,be)*(k(al)-p(al))); 
id vert(+3,-3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =    +g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))    
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,+3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =    -g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))     
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,n?{2,-2},+3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
id vert(+3,n?{2,-2},-3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
*
id vert(+3,-6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*pm(3)*d_(mu,be);
id vert(+3,n?{1,-1},-6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*pm(3)*d_(mu,al);
id vert(-3,+6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*pm(3)*d_(mu,be);
id vert(-3,n?{1,-1},+6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*pm(3)*d_(mu,al);
*
id vert(+3,-6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(+3,n?{2,-2},-6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(-3,+6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(-3,n?{2,-2},+6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,al);
*
id vert(+3,-3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =     -g*pm(3)*d_(mu,al);
id vert(-3,+3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =     -g*pm(3)*d_(mu,al);
id vert(+3,n?{4,-4},-3,mu?,al?,be?,p?,q?,k?) =     -g*pm(3)*d_(mu,be);
id vert(-3,n?{4,-4},+3,mu?,al?,be?,p?,q?,k?) =     -g*pm(3)*d_(mu,be);
*
id vert(+3,-6,n?{5,-5},mu?,al?,be?,p?,q?,k?) =  1/2*g*(q(mu)-k(mu));
id vert(-3,+6,n?{5,-5},mu?,al?,be?,p?,q?,k?) =  1/2*g*(k(mu)-q(mu));
id vert(+3,n?{5,-5},-6,mu?,al?,be?,p?,q?,k?) =  1/2*g*(k(mu)-q(mu));
id vert(-3,n?{5,-5},+6,mu?,al?,be?,p?,q?,k?) =  1/2*g*(q(mu)-k(mu));
*
id vert(+3,-6,n?{4,-4},mu?,al?,be?,p?,q?,k?) = i_/2*g*(q(mu)-k(mu));
id vert(-3,+6,n?{4,-4},mu?,al?,be?,p?,q?,k?) = i_/2*g*(q(mu)-k(mu));
id vert(+3,n?{4,-4},-6,mu?,al?,be?,p?,q?,k?) = i_/2*g*(k(mu)-q(mu));
id vert(-3,n?{4,-4},+6,mu?,al?,be?,p?,q?,k?) = i_/2*g*(k(mu)-q(mu));
*
id vert(+3,-9,+8,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(2)*q(mu);
id vert(+3,+8,-9,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(2)*k(mu);
id vert(-3,-8,+9,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(3)*q(mu);
id vert(-3,+9,-8,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(3)*k(mu);
id vert(+3,-10,+8,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(1)*q(mu);
id vert(+3,+8,-10,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(1)*k(mu);
id vert(-3,-8,+10,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(3)*q(mu);
id vert(-3,+10,-8,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(3)*k(mu);
*
id vert(+3,-7,+9,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(3)*q(mu);
id vert(+3,+9,-7,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(3)*k(mu);
id vert(-3,-9,+7,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(2)*q(mu);
id vert(-3,+7,-9,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(2)*k(mu);
*
id vert(+3,-7,+10,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(3)*q(mu);
id vert(+3,+10,-7,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(3)*k(mu);
id vert(-3,-10,+7,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(1)*q(mu);
id vert(-3,+7,-10,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(1)*k(mu);
*
id vert(n1?{4,-4},n2?{2,-2},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = -g*pm(3)/ctw^2*d_(al,be);
id vert(n1?{4,-4},n2?{2,-2},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*(k(al)-p(al));
id vert(n1?{4,-4},n2?{5,-5},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*(q(be)-p(be));
id vert(n1?{4,-4},n2?{4,-4},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = -3/2*g*pm(4)*pm(4)/pm(3);
id vert(n1?{4,-4},n2?{5,-5},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n1?{4,-4},+3,-3,mu?,al?,be?,p?,q?,k?)=                 -g*pm(3)*d_(al,be);
id vert(n1?{4,-4},-3,+3,mu?,al?,be?,p?,q?,k?)=                 -g*pm(3)*d_(al,be);
id vert(n?{4,-4},+3,-6,mu?,al?,be?,p?,q?,k?) =             i_/2*g*(k(al)-p(al));
id vert(n?{4,-4},-3,+6,mu?,al?,be?,p?,q?,k?) =             i_/2*g*(k(al)-p(al));
id vert(n?{4,-4},+6,-3,mu?,al?,be?,p?,q?,k?) =             i_/2*g*(q(be)-p(be));
id vert(n?{4,-4},-6,+3,mu?,al?,be?,p?,q?,k?) =             i_/2*g*(q(be)-p(be));
id vert(n?{4,-4},+6,-6,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{4,-4},-6,+6,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{4,-4},+7,-7,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},-7,+7,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},+8,-8,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},-8,+8,mu?,al?,be?,p?,q?,k?) =             -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},+9,-9,mu?,al?,be?,p?,q?,k?) =             -1/2*g/ctw/ctw*pm(3)*rxi(2);
id vert(n?{4,-4},-9,+9,mu?,al?,be?,p?,q?,k?) =             -1/2*g/ctw/ctw*pm(3)*rxi(2);
*
id vert(n1?{5,-5},n2?{2,-2},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*(p(al)-k(al));
id vert(n1?{5,-5},n2?{4,-4},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*(p(be)-q(be));
id vert(n1?{5,-5},n2?{5,-5},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n1?{5,-5},n2?{4,-4},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{5,-5},+3,-6,mu?,al?,be?,p?,q?,k?) =              1/2*g*(k(al)-p(al));
id vert(n?{5,-5},-3,+6,mu?,al?,be?,p?,q?,k?) =              1/2*g*(p(al)-k(al));
id vert(n?{5,-5},+6,-3,mu?,al?,be?,p?,q?,k?) =              1/2*g*(p(be)-q(be));
id vert(n?{5,-5},-6,+3,mu?,al?,be?,p?,q?,k?) =              1/2*g*(q(be)-p(be));
id vert(n?{5,-5},+7,-7,mu?,al?,be?,p?,q?,k?) =             i_/2*g*pm(3)*rxi(3);
id vert(n?{5,-5},-7,+7,mu?,al?,be?,p?,q?,k?) =             i_/2*g*pm(3)*rxi(3);
id vert(n?{5,-5},+8,-8,mu?,al?,be?,p?,q?,k?) =            -i_/2*g*pm(3)*rxi(3);
id vert(n?{5,-5},-8,+8,mu?,al?,be?,p?,q?,k?) =            -i_/2*g*pm(3)*rxi(3);
*
id vert(-6,+3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*pm(3)*d_(al,be);
id vert(-6,n?{1,-1},+3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*pm(3)*d_(be,al);
id vert(+6,-3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*pm(3)*d_(al,be);
id vert(+6,n?{1,-1},-3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*pm(3)*d_(be,al);
*
id vert(-6,+3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(al,be);
id vert(-6,n?{2,-2},+3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(be,al);
id vert(+6,-3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(al,be);
id vert(+6,n?{2,-2},-3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(be,al);
*
id vert(+6,-6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =     +g*stw*(p(be)-q(be));
id vert(-6,+6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =     +g*stw*(q(be)-p(be));
id vert(+6,n?{1,-1},-6,mu?,al?,be?,p?,q?,k?) =     +g*stw*(p(al)-k(al));
id vert(-6,n?{1,-1},+6,mu?,al?,be?,p?,q?,k?) =     +g*stw*(k(al)-p(al));
id vert(+6,-6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =     +g*(2*ctw-1/ctw)/2*(p(be)-q(be));
id vert(-6,+6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =     +g*(2*ctw-1/ctw)/2*(q(be)-p(be));
id vert(+6,n?{2,-2},-6,mu?,al?,be?,p?,q?,k?) =     +g*(2*ctw-1/ctw)/2*(p(al)-k(al));
id vert(-6,n?{2,-2},+6,mu?,al?,be?,p?,q?,k?) =     +g*(2*ctw-1/ctw)/2*(k(al)-p(al));
id vert(+6,-3,n?{5,-5},mu?,al?,be?,p?,q?,k?) =  1/2*g*(k(al)-p(al));
id vert(-6,+3,n?{5,-5},mu?,al?,be?,p?,q?,k?) =  1/2*g*(p(al)-k(al));
id vert(+6,n?{5,-5},-3,mu?,al?,be?,p?,q?,k?) =  1/2*g*(q(be)-p(be));
id vert(-6,n?{5,-5},+3,mu?,al?,be?,p?,q?,k?) =  1/2*g*(p(be)-q(be));
id vert(+6,-3,n?{4,-4},mu?,al?,be?,p?,q?,k?) = i_/2*g*(p(al)-k(al));
id vert(-6,+3,n?{4,-4},mu?,al?,be?,p?,q?,k?) = i_/2*g*(p(al)-k(al));
id vert(+6,n?{4,-4},-3,mu?,al?,be?,p?,q?,k?) = i_/2*g*(p(be)-q(be));
id vert(-6,n?{4,-4},+3,mu?,al?,be?,p?,q?,k?) = i_/2*g*(p(be)-q(be));
id vert(+6,-6,n?{4,-4},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(-6,+6,n?{4,-4},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(+6,n?{4,-4},-6,mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(-6,n?{4,-4},+6,mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
*
id vert(+6,-9,+8,mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*rxi(2)*pm(3);
id vert(+6,+8,-9,mu?,al?,be?,p?,q?,k?) = i_/2*g/ctw*rxi(2)*pm(3);
id vert(-6,-8,+9,mu?,al?,be?,p?,q?,k?) =   i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
id vert(-6,+9,-8,mu?,al?,be?,p?,q?,k?) =   i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
*
id vert(-6,-9,+7,mu?,al?,be?,p?,q?,k?) =-i_/2*g/ctw*rxi(2)*pm(3);
id vert(-6,+7,-9,mu?,al?,be?,p?,q?,k?) =-i_/2*g/ctw*rxi(2)*pm(3);
id vert(+6,-7,+9,mu?,al?,be?,p?,q?,k?) =  -i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
id vert(+6,+9,-7,mu?,al?,be?,p?,q?,k?) =  -i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
*
id vert(+6,-7,+10,mu?,al?,be?,p?,q?,k?)=  -i_*g*stw*rxi(3)*pm(3);
id vert(+6,+10,-7,mu?,al?,be?,p?,q?,k?)=  -i_*g*stw*rxi(3)*pm(3);
id vert(-6,-8,+10,mu?,al?,be?,p?,q?,k?)=   i_*g*stw*rxi(3)*pm(3);
id vert(-6,+10,-8,mu?,al?,be?,p?,q?,k?)=   i_*g*ctw*rxi(3)*pm(3);
*
*-             ----------boson_fermion_vertices-----------
*
* first generation
* ----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+13,-13,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-13,+13,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+14,-14,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-14,+14,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(-3,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(13)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(13)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(14)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(14)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii))*d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+11,-11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-11,+11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+12,-12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-12,+12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+13,-13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-13,+13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+14,-14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-14,+14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{2,-2},+11,-11,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},-11,+11,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},+12,-12,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},-12,+12,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},+13,-13,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},-13,+13,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},+14,-14,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
id vert(n?{2,-2},-14,+14,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
*
id vert(-3,+11,-12,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-12,+11,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+12,-11,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-11,+12,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+13,-14,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-14,+13,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+14,-13,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-13,+14,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
id vert(n?{4,-4},+11,-11,mu?,ii?) = -g/4*pm(11)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-11,+11,mu?,ii?) = -g/4*pm(11)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+12,-12,mu?,ii?) = -g/4*pm(12)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-12,+12,mu?,ii?) = -g/4*pm(12)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+13,-13,mu?,ii?) = -g/4*pm(13)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-13,+13,mu?,ii?) = -g/4*pm(13)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+14,-14,mu?,ii?) = -g/4*pm(14)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-14,+14,mu?,ii?) = -g/4*pm(14)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{5,-5},+11,-11,mu?,ii?) = i_*g*(vpa(11)-vma(11))/4/ctw*pm(11)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-11,+11,mu?,ii?) = i_*g*(vpa(11)-vma(11))/4/ctw*pm(11)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+12,-12,mu?,ii?) = i_*g*(vpa(12)-vma(12))/4/ctw*pm(12)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-12,+12,mu?,ii?) = i_*g*(vpa(12)-vma(12))/4/ctw*pm(12)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+13,-13,mu?,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-13,+13,mu?,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+14,-14,mu?,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-14,+14,mu?,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii));
*
id vert(-6,+11,-12,mu?,ii?) = i_*g/2/sr2*(pm(12)*pm(3)^-1*gd6(ii)-pm(11)*pm(3)^-1*gd7(ii));
id vert(-6,-12,+11,mu?,ii?) = i_*g/2/sr2*(pm(12)*pm(3)^-1*gd6(ii)-pm(11)*pm(3)^-1*gd7(ii));
id vert(+6,+12,-11,mu?,ii?) = i_*g/2/sr2*(pm(11)*pm(3)^-1*gd6(ii)-pm(12)*pm(3)^-1*gd7(ii));
id vert(+6,-11,+12,mu?,ii?) = i_*g/2/sr2*(pm(11)*pm(3)^-1*gd6(ii)-pm(12)*pm(3)^-1*gd7(ii));
id vert(-6,+13,-14,mu?,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii));
id vert(-6,-14,+13,mu?,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii));
id vert(+6,+14,-13,mu?,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii));
id vert(+6,-13,+14,mu?,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii));
*
* second generation
* -----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+17,-17,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-17,+17,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+18,-18,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-18,+18,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(-3,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(17)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(17)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(18)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(18)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii))*d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+15,-15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-15,+15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+16,-16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-16,+16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+17,-17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-17,+17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+18,-18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-18,+18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{2,-2},+15,-15,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},-15,+15,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},+16,-16,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},-16,+16,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},+17,-17,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},-17,+17,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},+18,-18,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
id vert(n?{2,-2},-18,+18,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
*
id vert(-3,+15,-16,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-16,+15,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+16,-15,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-15,+16,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+17,-18,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-18,+17,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+18,-17,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-17,+18,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
id vert(n?{4,-4},+15,-15,mu?,ii?) = -g/4*pm(15)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-15,+15,mu?,ii?) = -g/4*pm(15)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+16,-16,mu?,ii?) = -g/4*pm(16)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-16,+16,mu?,ii?) = -g/4*pm(16)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+17,-17,mu?,ii?) = -g/4*pm(17)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-17,+17,mu?,ii?) = -g/4*pm(17)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+18,-18,mu?,ii?) = -g/4*pm(18)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-18,+18,mu?,ii?) = -g/4*pm(18)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{5,-5},+15,-15,mu?,ii?) = i_*g*(vpa(15)-vma(15))/4/ctw*pm(15)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-15,+15,mu?,ii?) = i_*g*(vpa(15)-vma(15))/4/ctw*pm(15)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+16,-16,mu?,ii?) = i_*g*(vpa(16)-vma(16))/4/ctw*pm(16)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-16,+16,mu?,ii?) = i_*g*(vpa(16)-vma(16))/4/ctw*pm(16)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+17,-17,mu?,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-17,+17,mu?,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+18,-18,mu?,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-18,+18,mu?,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii));
*
id vert(-6,+15,-16,mu?,ii?) = i_*g/2/sr2*(pm(16)*pm(3)^-1*gd6(ii)-pm(15)*pm(3)^-1*gd7(ii));
id vert(-6,-16,+15,mu?,ii?) = i_*g/2/sr2*(pm(16)*pm(3)^-1*gd6(ii)-pm(15)*pm(3)^-1*gd7(ii));
id vert(+6,+16,-15,mu?,ii?) = i_*g/2/sr2*(pm(15)*pm(3)^-1*gd6(ii)-pm(16)*pm(3)^-1*gd7(ii));
id vert(+6,-15,+16,mu?,ii?) = i_*g/2/sr2*(pm(15)*pm(3)^-1*gd6(ii)-pm(16)*pm(3)^-1*gd7(ii));
id vert(-6,+17,-18,mu?,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii));
id vert(-6,-18,+17,mu?,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii));
id vert(+6,+18,-17,mu?,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii));
id vert(+6,-17,+18,mu?,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii));
*
* third generation
* ----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(n1)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(n2)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4/ctw*gd(ii,mu)*(vpa(n1)*gd6(ii)+vma(n1)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4/ctw*gd(ii,mu)*(vpa(n2)*gd6(ii)+vma(n2)*gd7(ii))*d(cl1,cl2);
id vert(-3,+21,-22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-22,+21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+22,-21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-21,+22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(n1)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{4,-4},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(n2)*pm(3)^-1*(gd6(ii)+gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(n1)-vma(n1))/4/ctw*pm(n1)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(n2)-vma(n2))/4/ctw*pm(n2)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+21,-22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-22,+21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+22,-21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-21,+22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
* SM vertices
id vert(n?{1,-1},+19,-19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-19,+19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+20,-20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-20,+20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+21,-21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-21,+21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},+22,-22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{1,-1},-22,+22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(gd6(ii)+gd7(ii));
id vert(n?{2,-2},+19,-19,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},-19,+19,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},+20,-20,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},-20,+20,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},+21,-21,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},-21,+21,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},+22,-22,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
id vert(n?{2,-2},-22,+22,mu?,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
*
id vert(-3,+19,-20,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-20,+19,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+20,-19,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-19,+20,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,+21,-22,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(-3,-22,+21,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,+22,-21,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
id vert(+3,-21,+22,mu?,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii);
*
id vert(n?{4,-4},+19,-19,mu?,ii?) = -g/4*pm(19)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-19,+19,mu?,ii?) = -g/4*pm(19)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+20,-20,mu?,ii?) = -g/4*pm(20)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-20,+20,mu?,ii?) = -g/4*pm(20)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+21,-21,mu?,ii?) = -g/4*pm(21)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-21,+21,mu?,ii?) = -g/4*pm(21)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},+22,-22,mu?,ii?) = -g/4*pm(22)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{4,-4},-22,+22,mu?,ii?) = -g/4*pm(22)*pm(3)^-1*(gd6(ii)+gd7(ii));
id vert(n?{5,-5},+19,-19,mu?,ii?) = i_*g*(vpa(19)-vma(19))/4/ctw*pm(19)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-19,+19,mu?,ii?) = i_*g*(vpa(19)-vma(19))/4/ctw*pm(19)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+20,-20,mu?,ii?) = i_*g*(vpa(20)-vma(20))/4/ctw*pm(20)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-20,+20,mu?,ii?) = i_*g*(vpa(20)-vma(20))/4/ctw*pm(20)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+21,-21,mu?,ii?) = i_*g*(vpa(21)-vma(21))/4/ctw*pm(21)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-21,+21,mu?,ii?) = i_*g*(vpa(21)-vma(21))/4/ctw*pm(21)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},+22,-22,mu?,ii?) = i_*g*(vpa(22)-vma(22))/4/ctw*pm(22)*pm(2)^-1*(gd6(ii)-gd7(ii));
id vert(n?{5,-5},-22,+22,mu?,ii?) = i_*g*(vpa(22)-vma(22))/4/ctw*pm(22)*pm(2)^-1*(gd6(ii)-gd7(ii));
*
id vert(-6,+19,-20,mu?,ii?) = i_*g/2/sr2*(pm(20)*pm(3)^-1*gd6(ii)-pm(19)*pm(3)^-1*gd7(ii));
id vert(-6,-20,+19,mu?,ii?) = i_*g/2/sr2*(pm(20)*pm(3)^-1*gd6(ii)-pm(19)*pm(3)^-1*gd7(ii));
id vert(+6,+20,-19,mu?,ii?) = i_*g/2/sr2*(pm(19)*pm(3)^-1*gd6(ii)-pm(20)*pm(3)^-1*gd7(ii));
id vert(+6,-19,+20,mu?,ii?) = i_*g/2/sr2*(pm(19)*pm(3)^-1*gd6(ii)-pm(20)*pm(3)^-1*gd7(ii));
id vert(-6,+21,-22,mu?,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii));
id vert(-6,-22,+21,mu?,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii));
id vert(+6,+22,-21,mu?,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii));
id vert(+6,-21,+22,mu?,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii));
*    
*           ---------boson_boson_four_line_vertices-----------
*
*- 3311: 12 Rules
id vert(n1?{1,-1},n2?{1,-1},+3,-3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n1?{1,-1},n2?{1,-1},-3,+3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(n1?{1,-1},+3,-3,n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n1?{1,-1},-3,+3,n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(+3,-3,n1?{1,-1},n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(-3,+3,n1?{1,-1},n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
*
id vert(n1?{1,-1},+3,n2?{1,-1},-3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(al,nu)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(n1?{1,-1},-3,n2?{1,-1},+3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(al,nu)
                                                                  -d_(mu,nu)*d_(al,be));
*
id vert(+3,n1?{1,-1},n2?{1,-1},-3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(-3,n1?{1,-1},n2?{1,-1},+3,mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(+3,n1?{1,-1},-3,n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(-3,n1?{1,-1},+3,n2?{1,-1},mu?,nu?,al?,be?) = -g^2*stw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
*- 3322: 12 Rules
id vert(n1?{2,-2},n2?{2,-2},+3,-3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n1?{2,-2},n2?{2,-2},-3,+3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(n1?{2,-2},+3,-3,n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n1?{2,-2},-3,+3,n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,be)*d_(al,nu)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(+3,-3,n1?{2,-2},n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(-3,+3,n1?{2,-2},n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
*-
id vert(n1?{2,-2},+3,n2?{2,-2},-3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(n1?{2,-2},-3,n2?{2,-2},+3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
*
id vert(+3,n1?{2,-2},n2?{2,-2},-3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(-3,n1?{2,-2},n2?{2,-2},+3,mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*
id vert(+3,n1?{2,-2},-3,n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(-3,n1?{2,-2},+3,n2?{2,-2},mu?,nu?,al?,be?) = -g^2*ctw^2*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
*- 3312: 24 Rules
id vert(n1?{1,-1},n2?{2,-2},+3,-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(n1?{1,-1},n2?{2,-2},-3,+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(n1?{2,-2},n2?{1,-1},+3,-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(n1?{2,-2},n2?{1,-1},-3,+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
*-
id vert(+3,-3,n1?{1,-1},n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(-3,+3,n1?{1,-1},n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(+3,-3,n1?{2,-2},n2?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
id vert(-3,+3,n1?{2,-2},n2?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al));
*-
id vert(n1?{1,-1},+3,n2?{2,-2},-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al) 
                                                                  -d_(mu,nu)*d_(al,be));
id vert(n1?{1,-1},-3,n2?{2,-2},+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(n2?{2,-2},+3,n1?{1,-1},-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al) 
                                                                  -d_(mu,nu)*d_(al,be));
id vert(n2?{2,-2},-3,n1?{1,-1},+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
*-
id vert(+3,n1?{1,-1},-3,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al) 
                                                                  -d_(mu,nu)*d_(al,be));
id vert(-3,n1?{1,-1},+3,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
id vert(+3,n2?{2,-2},-3,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al) 
                                                                  -d_(mu,nu)*d_(al,be));
id vert(-3,n2?{2,-2},+3,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,al)*d_(nu,be)
                                                                  -d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be));
*-
id vert(n1?{1,-1},+3,-3,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be) 
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n1?{1,-1},-3,+3,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n2?{2,-2},+3,-3,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be) 
                                                                  -d_(mu,al)*d_(nu,be));
id vert(n2?{2,-2},-3,+3,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*-
id vert(+3,n1?{1,-1},n2?{2,-2},-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be) 
                                                                  -d_(mu,al)*d_(nu,be));
id vert(-3,n1?{1,-1},n2?{2,-2},+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
id vert(+3,n2?{2,-2},n1?{1,-1},-3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be) 
                                                                  -d_(mu,al)*d_(nu,be));
id vert(-3,n2?{2,-2},n1?{1,-1},+3,mu?,nu?,al?,be?)=-g^2*stw*ctw*(2*d_(mu,be)*d_(nu,al)
                                                                  -d_(mu,nu)*d_(al,be)
                                                                  -d_(mu,al)*d_(nu,be));
*-
*-6611: 12 Rules
id vert(n1?{1,-1},n2?{1,-1},+6,-6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,nu);
id vert(n1?{1,-1},n2?{1,-1},-6,+6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,nu);
id vert(n1?{1,-1},+6,-6,n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,be);
id vert(n1?{1,-1},-6,+6,n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,be);
id vert(+6,-6,n1?{1,-1},n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(al,be);
id vert(-6,+6,n1?{1,-1},n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(al,be);
id vert(n1?{1,-1},+6,n2?{1,-1},-6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,al);
id vert(n1?{1,-1},-6,n2?{1,-1},+6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(mu,al);
id vert(+6,n1?{1,-1},-6,n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(nu,be);
id vert(-6,n1?{1,-1},+6,n2?{1,-1},mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(nu,be);
id vert(+6,n1?{1,-1},n2?{1,-1},-6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(nu,al);
id vert(-6,n1?{1,-1},n2?{1,-1},+6,mu?,nu?,al?,be?) = -2*g^2*stw^2*d_(nu,al);
*-6622: 12 Rules
id vert(n1?{2,-2},n2?{2,-2},+6,-6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,nu);
id vert(n1?{2,-2},n2?{2,-2},-6,+6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,nu);
id vert(n1?{2,-2},+6,-6,n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,be);
id vert(n1?{2,-2},-6,+6,n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,be);
id vert(+6,-6,n1?{2,-2},n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(al,be);
id vert(-6,+6,n1?{2,-2},n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(al,be);
id vert(n1?{2,-2},+6,n2?{2,-2},-6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,al);
id vert(n1?{2,-2},-6,n2?{2,-2},+6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(mu,al);
id vert(+6,n1?{2,-2},-6,n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(nu,be);
id vert(-6,n1?{2,-2},+6,n2?{2,-2},mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(nu,be);
id vert(+6,n1?{2,-2},n2?{2,-2},-6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(nu,al);
id vert(-6,n1?{2,-2},n2?{2,-2},+6,mu?,nu?,al?,be?)=-1/2*g^2/ctw^2*(2*ctw^2-1)^2*d_(nu,al);
*-6612: 24 Rules
id vert(n1?{1,-1},n2?{2,-2},+6,-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,nu);
id vert(n1?{1,-1},n2?{2,-2},-6,+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,nu);
id vert(n1?{1,-1},+6,n2?{2,-2},-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,al);
id vert(n1?{1,-1},-6,n2?{2,-2},+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,al);
id vert(n1?{1,-1},+6,-6,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,be);
id vert(n1?{1,-1},-6,+6,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,be);
id vert(+6,n1?{1,-1},-6,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,be);
id vert(-6,n1?{1,-1},+6,n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,be);
id vert(+6,n1?{1,-1},n2?{2,-2},-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,al);
id vert(-6,n1?{1,-1},n2?{2,-2},+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,al);
id vert(+6,-6,n1?{1,-1},n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(al,be);
id vert(-6,+6,n1?{1,-1},n2?{2,-2},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(al,be);
id vert(n2?{2,-2},n1?{1,-1},+6,-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,nu);
id vert(n2?{2,-2},n1?{1,-1},-6,+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,nu);
id vert(n2?{2,-2},+6,n1?{1,-1},-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,al);
id vert(n2?{2,-2},-6,n1?{1,-1},+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,al);
id vert(n2?{2,-2},+6,-6,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,be);
id vert(n2?{2,-2},-6,+6,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(mu,be);
id vert(+6,n2?{2,-2},-6,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,be);
id vert(-6,n2?{2,-2},+6,n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,be);
id vert(+6,n2?{2,-2},n1?{1,-1},-6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,al);
id vert(-6,n2?{2,-2},n1?{1,-1},+6,mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(nu,al);
id vert(+6,-6,n2?{2,-2},n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(al,be);
id vert(-6,+6,n2?{2,-2},n1?{1,-1},mu?,nu?,al?,be?)=-g^2*stw/ctw*(2*ctw^2-1)*d_(al,be);
* t2 4 1 63/36 
id vert(n1?{4,-4},n2?{1,-1},+6,-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},n2?{1,-1},-6,+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},n2?{1,-1},+3,-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(n1?{4,-4},n2?{1,-1},-3,+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
* t2 4 63 1 36                   
id vert(n1?{4,-4},+6,n2?{1,-1},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},-6,n2?{1,-1},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},+3,n2?{1,-1},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(n1?{4,-4},-3,n2?{1,-1},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
* t2 4 63/36 1
id vert(n1?{4,-4},+6,-3,n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},-6,+3,n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},+3,-6,n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},-3,+6,n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
* t2 63/36 41
id vert(+6,-3,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
id vert(-6,+3,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(+3,-6,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,be);
id vert(-3,+6,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,be);
* t2 63/36 14
id vert(+6,-3,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
id vert(-6,+3,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(+3,-6,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(-3,+6,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
* t2 1 63/36 4  New 14 July 2004  !!! to be continued
id vert(n1?{1,-1},+6,-3,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},-6,+3,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},+3,-6,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,nu);
id vert(n1?{1,-1},-3,+6,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,nu);
* t2 63/36 51
id vert(+6,-3,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(-6,+3,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(+3,-6,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(-3,+6,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
* t2 63/36 15
id vert(+6,-3,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(-6,+3,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(+3,-6,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
id vert(-3,+6,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
* t2 51 63/36
id vert(n1?{5,-5},n2?{1,-1},+6,-3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(n1?{5,-5},n2?{1,-1},-6,+3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(n1?{5,-5},n2?{1,-1},+3,-6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(n1?{5,-5},n2?{1,-1},-3,+6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
* t2 15 63/36
id vert(n1?{1,-1},n2?{5,-5},+6,-3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(n1?{1,-1},n2?{5,-5},-6,+3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(n1?{1,-1},n2?{5,-5},+3,-6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},n2?{5,-5},-3,+6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
* t2 52 63/36
id vert(n1?{5,-5},n2?{2,-2},+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{5,-5},n2?{2,-2},-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{5,-5},n2?{2,-2},+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,al);
id vert(n1?{5,-5},n2?{2,-2},-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,al);
* t2 25 63/36
id vert(n1?{2,-2},n2?{5,-5},+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},n2?{5,-5},-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,al);
* t7 hgg
id vert(-3,n1?{4,-4},n2?{1,-1},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
id vert(+3,n1?{4,-4},n2?{1,-1},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(-6,n1?{4,-4},n2?{1,-1},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(+6,n1?{4,-4},n2?{1,-1},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
*
id vert(-3,n2?{1,-1},n1?{4,-4},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,nu);
id vert(+3,n2?{1,-1},n1?{4,-4},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,nu);
id vert(-6,n2?{1,-1},n1?{4,-4},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(+6,n2?{1,-1},n1?{4,-4},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
*-2634 dobavlyat'....
id vert(n1?{2,-2},+6,-3,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},-6,+3,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},+3,-6,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,nu);
id vert(n1?{2,-2},-3,+6,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,nu);
*-4263
id vert(n1?{4,-4},n2?{2,-2},+6,-3,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{4,-4},n2?{2,-2},-6,+3,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{4,-4},n2?{2,-2},+3,-6,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,al);
id vert(n1?{4,-4},n2?{2,-2},-3,+6,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,al);
*
id vert(n1?{4,-4},+6,n2?{2,-2},-3,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(al,be);
id vert(n1?{4,-4},-6,n2?{2,-2},+3,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(al,be);
id vert(n1?{4,-4},+3,n2?{2,-2},-6,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,al);
id vert(n1?{4,-4},-3,n2?{2,-2},+6,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,al);
*
id vert(n1?{4,-4},+6,-3,n2?{2,-2},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(al,be);
id vert(n1?{4,-4},-6,+3,n2?{2,-2},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(al,be);
id vert(n1?{4,-4},+3,-6,n2?{2,-2},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{4,-4},-3,+6,n2?{2,-2},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,be);
*
id vert(+6,-3,n1?{4,-4},n2?{2,-2},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(-6,+3,n1?{4,-4},n2?{2,-2},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(+3,-6,n1?{4,-4},n2?{2,-2},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,be);
id vert(-3,+6,n1?{4,-4},n2?{2,-2},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,be);
* 
id vert(-3,n1?{4,-4},n2?{2,-2},+6,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(+3,n1?{4,-4},n2?{2,-2},-6,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(-6,n1?{4,-4},n2?{2,-2},+3,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(al,be);
id vert(+6,n1?{4,-4},n2?{2,-2},-3,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(al,be);
*			          		 
id vert(-3,n2?{2,-2},n1?{4,-4},+6,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,nu);
id vert(+3,n2?{2,-2},n1?{4,-4},-6,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,nu);
id vert(-6,n2?{2,-2},n1?{4,-4},+3,mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(nu,be);
id vert(+6,n2?{2,-2},n1?{4,-4},-3,mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(nu,be);
*--
id vert(n1?{2,-2},n2?{2,-2},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,nu);
id vert(n1?{2,-2},n2?{4,-4},n3?{4,-4},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,be);
id vert(n1?{2,-2},n2?{4,-4},n3?{2,-2},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,al);
id vert(n1?{4,-4},n2?{2,-2},n3?{2,-2},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,al);
id vert(n1?{4,-4},n2?{2,-2},n3?{4,-4},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,be);
id vert(n1?{4,-4},n2?{4,-4},n3?{2,-2},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw/ctw*d_(al,be);
*-
id vert(n1?{2,-2},n2?{2,-2},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,nu);
id vert(n1?{2,-2},n2?{5,-5},n3?{5,-5},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},n3?{2,-2},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,al);
id vert(n1?{5,-5},n2?{2,-2},n3?{2,-2},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,al);
id vert(n1?{5,-5},n2?{2,-2},n3?{5,-5},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,be);
id vert(n1?{5,-5},n2?{5,-5},n3?{2,-2},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(al,be);
*-
id vert(+3,-3,+3,-3,mu?,nu?,al?,be?)= +g^2*(2*d_(mu,al)*d_(nu,be)  
				             -d_(mu,nu)*d_(al,be)  
				             -d_(mu,be)*d_(nu,al));
id vert(-3,+3,-3,+3,mu?,nu?,al?,be?)= +g^2*(2*d_(mu,al)*d_(nu,be)  
				             -d_(mu,nu)*d_(al,be)  
				             -d_(mu,be)*d_(nu,al));
id vert(-3,-3,+3,+3,mu?,nu?,al?,be?)= +g^2*(2*d_(mu,nu)*d_(al,be)
                                             -d_(mu,al)*d_(nu,be)
                                             -d_(mu,be)*d_(nu,al));
id vert(-3,+3,+3,-3,mu?,nu?,al?,be?)= +g^2*(2*d_(mu,be)*d_(nu,al)
                                             -d_(mu,nu)*d_(al,be)
                                             -d_(mu,al)*d_(nu,be));
id vert(+3,-3,-3,+3,mu?,nu?,al?,be?)= +g^2*(2*d_(nu,al)*d_(mu,be)
                                             -d_(nu,mu)*d_(al,be)
                                             -d_(nu,be)*d_(al,mu));
id vert(+3,+3,-3,-3,mu?,nu?,al?,be?)= +g^2*(2*d_(nu,mu)*d_(al,be)
                                             -d_(nu,al)*d_(mu,be)
                                             -d_(nu,be)*d_(al,mu));
*-
id vert(-3,+3,n1?{4,-4},n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,nu);
id vert(+3,-3,n1?{4,-4},n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,mu);
id vert(n1?{4,-4},n2?{4,-4},+3,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(n1?{4,-4},n2?{4,-4},-3,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,al);
id vert(-3,n1?{4,-4},n2?{4,-4},+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(-3,n1?{4,-4},+3,n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
id vert(n1?{4,-4},-3,n2?{4,-4},+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(n1?{4,-4},-3,+3,n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(+3,n1?{4,-4},n2?{4,-4},-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,mu);
id vert(+3,n1?{4,-4},-3,n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(al,mu);
id vert(n1?{4,-4},+3,n2?{4,-4},-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,nu);
id vert(n1?{4,-4},+3,-3,n2?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
*-
id vert(-3,+3,n1?{5,-5},n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,nu);
id vert(+3,-3,n1?{5,-5},n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,mu);
id vert(n1?{5,-5},n2?{5,-5},+3,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(n1?{5,-5},n2?{5,-5},-3,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,al);
id vert(-3,n1?{5,-5},n2?{5,-5},+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(-3,n1?{5,-5},+3,n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
id vert(n1?{5,-5},-3,n2?{5,-5},+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(n1?{5,-5},-3,+3,n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(+3,n1?{5,-5},n2?{5,-5},-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,mu);
id vert(+3,n1?{5,-5},-3,n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(al,mu);
id vert(n1?{5,-5},+3,n2?{5,-5},-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(be,nu);
id vert(n1?{5,-5},+3,-3,n2?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
*-
id vert(+3,-3,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,nu);
id vert(-3,+3,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,mu);
id vert(+3,-3,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,nu);
id vert(-3,+3,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,mu);
id vert(+6,-6,+3,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(+6,-6,-3,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(-6,+6,+3,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(-6,+6,-3,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(al,be);
id vert(+3,-6,+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(+3,-6,-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
id vert(-6,+3,+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(-6,+3,-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(+6,-3,+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(+6,-3,-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(-3,+6,+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
id vert(-3,+6,-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(+3,+6,-3,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
id vert(+3,+6,-6,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(+6,+3,-3,-6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(+6,+3,-6,-3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(-6,-3,+6,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,be);
id vert(-6,-3,+3,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(nu,al);
id vert(-3,-6,+6,+3,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,be);
id vert(-3,-6,+3,+6,mu?,nu?,al?,be?)= -1/2*g^2*d_(mu,al);
*-
id vert(n1?{4,-4},n2?{4,-4},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)=-3/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
*-
id vert(n1?{4,-4},n2?{4,-4},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},n2?{5,-5},n3?{5,-5},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},n2?{5,-5},n3?{4,-4},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{4,-4},n3?{5,-5},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{4,-4},n3?{4,-4},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{5,-5},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
*-
id vert(n1?{5,-5},n2?{5,-5},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)=-3/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
*-
id vert(n1?{4,-4},n2?{4,-4},+6,-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},n2?{4,-4},-6,+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},-6,+6,n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},-6,n2?{4,-4},+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,n1?{4,-4},+6,n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,n1?{4,-4},n2?{4,-4},+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,-6,n1?{4,-4},n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,n1?{4,-4},n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,n1?{4,-4},n2?{4,-4},-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,n1?{4,-4},-6,n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},+6,n2?{4,-4},-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},+6,-6,n2?{4,-4},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
*-
id vert(n1?{5,-5},n2?{5,-5},+6,-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{5,-5},-6,+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},-6,+6,n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},-6,n2?{5,-5},+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,n1?{5,-5},+6,n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,n1?{5,-5},n2?{5,-5},+6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,-6,n1?{5,-5},n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,n1?{5,-5},n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,n1?{5,-5},n2?{5,-5},-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,n1?{5,-5},-6,n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},+6,n2?{5,-5},-6,mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},+6,-6,n2?{5,-5},mu?,nu?,al?,be?)= -1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
*-
id vert(+6,-6,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,-6,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,-6,+6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,+6,-6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(?a)=0;
id sr2^-2 = 1/2;
*
*
* ###################################################
* #            P R O P A G A T O R S                #
* ###################################################
*
*   
* ========== SM Bosons ==========
* ===============================
*
id pr(n?{1,-1},mu?,nu?,p?)  = den(1,0,p)*(d_(mu,nu) + (rxi(1)^2-1)*p(mu)*p(nu)*den(1,0,p));
id pr(n?{2,-2},mu?,nu?,p?)  = den(1,pm(2),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(2)^-1*pm(2)^-1)
                              - p(mu)*p(nu)*pm(2)^-1*pm(2)^-1*den(rxi(2),pm(2),p);
id pr(n?{3,-3},mu?,nu?,p?)  = den(1,pm(3),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(3)^-1*pm(3)^-1)
                              - p(mu)*p(nu)*pm(3)^-1*pm(3)^-1*den(rxi(3),pm(3),p);
id pr(n?{4,-4},mu?,nu?,p?)  = den(1,pm(4),p);
id pr(n?{5,-5},mu?,nu?,p?)  = den(rxi(2),pm(2),p);
id pr(n?{6,-6},mu?,nu?,p?)  = den(rxi(3),pm(3),p);
id pr(n?{7,-7},mu?,nu?,p?)  = rxi(3)*den(rxi(3),pm(3),p);
id pr(n?{8,-8},mu?,nu?,p?)  = rxi(3)*den(rxi(3),pm(3),p);
id pr(n?{9,-9},mu?,nu?,p?)  = rxi(2)*den(rxi(2),pm(2),p);
id pr(n?{10,-10},mu?,nu?,p?)= rxi(1)*den(1,0,p);
*
id pr(n?{4,-4},p?)  = den(1,pm(4),p);
id pr(n?{5,-5},p?)  = den(rxi(2),pm(2),p);
id pr(n?{6,-6},p?)  = den(rxi(3),pm(3),p);
id pr(n?{7,-7},p?)  = rxi(3)*den(rxi(3),pm(3),p);
id pr(n?{8,-8},p?)  = rxi(3)*den(rxi(3),pm(3),p);
id pr(n?{9,-9},p?)  = rxi(2)*den(rxi(2),pm(2),p);
id pr(n?{10,-10},p?)= rxi(1)*den(1,0,p);
*
*
* ======= QCD Gluon/Ghost =======
* ===============================
*
id pr(23,mu?,nu?,gi1?sAgi,gi2?sAgi,p?)= d(gi1,gi2)*
                                        den(1,0,p)*(d_(mu,nu) + (rxi(23)^2 - 1)*p(mu)*p(nu)*den(1,0,p));
id pr(24,mu?,nu?,gi1?sAgi,gi2?sAgi,p?)= d(gi1,gi2)*den(1,0,p);
*
*
* == Fermions - Quarks/leptons ==
* ===============================
*
id pr(n?sQI,?a,cl1?scl,cl2?scl,i?) = d(cl1,cl2)*den(1,pm(abs_(n)),?a)*pr(sf,n,?a,i);
id pr(n?sFI,p?,i?) = den(1,pm(abs_(n)),p)*pr(sf,n,p,i);

splitarg pr 3;
repeat id pr(sf,n?sFI,p?,?a,i?) = -i_*gd(i,p)+pr(sf,n,?a,i); 
id pr(sf,n?sFI,i?) = pm(abs_(n));
id gd(i?,-p?) = -gd(i,p);
*
*
*
* ###################################################
* #           COMBINATORICAL FACTORS                #
* ###################################################
*
*
* ==== Three Point Vertices ====
* ==============================
*
id comb3(2,2, 4 ,2) = 0;
id comb3(2,2,-4,-2) = 0;
id comb3(2,2,+4,-2) = 0;
id comb3(2,2,-4,+2) = 0;
*-
id comb3(2,2, 4 ,5) = 0;
id comb3(2,2,-4,-5) = 0;
id comb3(2,2,-4,+5) = 0;
id comb3(2,2,+4,-5) = 0;
*-
id comb3(4,4, 2 ,2) = 1/2;
id comb3(4,4, 4, 4) = 1/2;
id comb3(4,4, 5, 5) = 1/2;
*-
id comb3(4,4,-2,-2) = 0;
id comb3(4,4,-2,+2) = 0;
id comb3(4,4,+2,-2) = 0;
*-
id comb3(4,4,-4,-4) = 0;
id comb3(4,4,+4,-4) = 0;
id comb3(4,4,-4,+4) = 0;
*-
id comb3(4,4,-5,-5) = 0;
id comb3(4,4,+5,-5) = 0;
id comb3(4,4,-5,+5) = 0;
*-
id comb3(4,4, 5, 2) = 0;
id comb3(4,4,-5,-2) = 0;
id comb3(4,4,+5,-2) = 0;
id comb3(4,4,-5,+2) = 0;
*-
id comb3(5,5, 4, 2) = 0;
id comb3(5,5,-4,-2) = 0;
id comb3(5,5,+4,-2) = 0;
id comb3(5,5,-4,+2) = 0;
*-
id comb3(5,5, 4, 5) = 0;
id comb3(5,5,-4,-5) = 0;
id comb3(5,5,-4,+5) = 0;
id comb3(5,5,+4,-5) = 0;
*- QCD -*
id comb3(23,23,23,23) = 1/2;
id comb3(23,23,24,24) = 1/2;
*-
id comb3(?a) = 1;
*
* ===== Four Point Vertices ====
* ==============================
*
id comb4( 2, 2,4) = +1/2;
id comb4( 2, 2,5) = +1/2;
id comb4( 2,-2,4) = +1/2;
id comb4( 2,-2,5) = +1/2;
id comb4(+3,-3,1) = +1/2;
id comb4(-3,+3,1) = +1/2;
id comb4(+3,-3,2) = +1/2;
id comb4(-3,+3,2) = +1/2;
id comb4(+3,-3,4) = +1/2;
id comb4(-3,+3,4) = +1/2;
id comb4(+3,-3,5) = +1/2;
id comb4(-3,+3,5) = +1/2;
id comb4( 4, 4,2) = +1/2;
id comb4( 4,-4,2) = +1/2;
id comb4( 4, 4,4) = +1/2;
id comb4( 4,-4,4) = +1/2;
id comb4( 4, 4,5) = +1/2;
id comb4( 4,-4,5) = +1/2;
id comb4( 5, 5,2) = +1/2;
id comb4( 5,-5,2) = +1/2;
id comb4( 5, 5,4) = +1/2;
id comb4( 5,-5,4) = +1/2;
id comb4( 5, 5,5) = +1/2;
id comb4( 5,-5,5) = +1/2;
id comb4(+6,-6,1) = +1/2;
id comb4(-6,+6,1) = +1/2;
id comb4(+6,-6,2) = +1/2;
id comb4(-6,+6,2) = +1/2;
id comb4(+6,-6,4) = +1/2;
id comb4(-6,+6,4) = +1/2;
id comb4(+6,-6,5) = +1/2;
id comb4(-6,+6,5) = +1/2;
*- QCD -*
id comb4(-23,+23, 23) = +1/2;
id comb4(+23,-23, 23) = +1/2;
*-
id comb4(?a) = 1;
*
* = tadpoles and for Higgs self energies =
* ========================================
*
id cft(2) = 1/2;
id cft(4) = 1/2;
id cft(5) = 1/2;
id cft(?a)= 1;
*
*
* ###################################################
* #            F I E L D   S I G N S                #
* ###################################################
*
#do i=7,22
 id sign(`i')= -1;
#enddo
 id sign(24) = -1;
 id sign(?a) =  1;
*
*
.sort

 repeat;
   id Tf(gi1?sAgi,?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = Tf(gi1,?a,?c,?d,?b);
   id d(?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = d(?a,?c,?d,?b);

   id Ta(?a,gi1?sAgi,?b)*d(?c,gi1?sAgi,?d) = Ta(?a,?c,?d,?b);
   id d(?a,cl1?sAgi,?b)*d(?c,cl1?sAgi,?d) = d(?a,?c,?d,?b);
 endrepeat;

#if {`I'} = 0
    #call a2b(gd6,gd5)
    #call a2b(gd7,gd5)
    #call a2b(g,e)
#endif
.sort
#endprocedure
*------------
