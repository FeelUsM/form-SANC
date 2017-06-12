#procedure FeynmanRules(I)
*----------------------------
* I=0: QED,QCD
* I=1: EW

.sort :Feynman-Start;
#include Misc.h
vector p,q,k; * в трехбозонных вершинах
Sym n1,n2;
index i;
Cfun	sf;	* for internal use
#include Declare.h
#call Globals()



* ###########################################################################################
* #                                                                                         #
* #                                FEYNMAN RULES FOR SM                                     #
* #                                                                                         #
* ###########################################################################################

*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             1 = gamma            2 = z           +/-3 = w^{+/-}    *
*             4 = h                5 = phi0        +/-6 = phi^{+/-}  *
*             7 = Xp               8 = Xm             9 = Yz         *
*             10= Ya               23= g              24= Xg         *
*====================================================================*
*             -------------boson_boson_vertices-----------
* --- 1 3 3 --- gamma W W ---
id vert(n?{1,-1},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(mu,al)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
id vert(n?{1,-1},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(mu,al)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
id vert(+3,-3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(-3,+3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(be,al)*(q(mu)-k(mu))  
                                                          +d_(be,mu)*(k(al)-p(al)));
id vert(+3,n?{1,-1},-3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))    
                                                          +d_(mu,be)*(k(al)-p(al))); 
id vert(-3,n?{1,-1},+3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(+d_(al,mu)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))      
                                                          +d_(mu,be)*(k(al)-p(al))); 
* --- 2 3 3 --- Z W W ---
id vert(n?{2,-2},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(+d_(mu,al)*(p(be)-q(be))
                                                          +d_(al,be)*(q(mu)-k(mu))
                                                          +d_(mu,be)*(k(al)-p(al)));
id vert(n?{2,-2},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(+d_(mu,al)*(p(be)-q(be))
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
* ----------------------------------
* --- 1 3 6 --- gamma W phi+- ---
id vert(n?{1,-1},-3,+6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(mu,al);
id vert(n?{1,-1},+6,-3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(mu,be);
id vert(n?{1,-1},+3,-6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(mu,al);
id vert(n?{1,-1},-6,+3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(mu,be);
id vert(+3,-6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(mu,be);
id vert(+3,n?{1,-1},-6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(mu,al);
id vert(-3,+6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(mu,be);
id vert(-3,n?{1,-1},+6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(mu,al);
id vert(-6,+3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(al,be);
id vert(-6,n?{1,-1},+3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw*      pm(3)*d_(be,al);
id vert(+6,-3,n?{1,-1},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(al,be);
id vert(+6,n?{1,-1},-3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw*      pm(3)*d_(be,al);
* --- 2 3 6 --- Z W phi+- ---
id vert(n?{2,-2},-3,+6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(n?{2,-2},+6,-3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(n?{2,-2},-6,+3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(n?{2,-2},+3,-6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(+3,-6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(+3,n?{2,-2},-6,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(-3,+6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,be);
id vert(-3,n?{2,-2},+6,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(mu,al);
id vert(-6,+3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(al,be);
id vert(-6,n?{2,-2},+3,mu?,al?,be?,p?,q?,k?) =  -i_*g*stw^2/ctw*pm(3)*d_(be,al);
id vert(+6,-3,n?{2,-2},mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(al,be);
id vert(+6,n?{2,-2},-3,mu?,al?,be?,p?,q?,k?) =  +i_*g*stw^2/ctw*pm(3)*d_(be,al);
* --- 2 2 4 --- Z Z h ---
id vert(n?{2,-2},n2?{2,-2},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = -g*pm(3)/ctw^2*d_(mu,al);
id vert(n?{2,-2},n2?{4,-4},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = -g*pm(3)/ctw^2*d_(mu,be);
id vert(n?{4,-4},n2?{2,-2},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = -g*pm(3)/ctw^2*d_(al,be);
* --- 3 3 4 --- W W h ---
id vert(+3,-3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(mu,al);
id vert(-3,+3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(mu,al);
id vert(+3,n?{4,-4},-3,mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(mu,be);
id vert(-3,n?{4,-4},+3,mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(mu,be);
id vert(n?{4,-4},+3,-3,mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(al,be);
id vert(n?{4,-4},-3,+3,mu?,al?,be?,p?,q?,k?) =               -g*pm(3)*      d_(al,be);
* ----------------------------------
* --- 4 4 4 --- h h h ---
id vert(n?{4,-4},n2?{4,-4},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = -3/2*g*pm(4)*pm(4)/pm(3);
* --- 4 5 5 --- h phi0 phi0 ---
id vert(n?{4,-4},n2?{5,-5},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{5,-5},n2?{5,-5},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{5,-5},n2?{4,-4},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = -1/2*g*pm(4)*pm(4)/pm(3);
* --- 4 6 6 --- 
id vert(n?{4,-4},+6,-6,mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
id vert(n?{4,-4},-6,+6,mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
id vert(+6,-6,n?{4,-4},mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
id vert(-6,+6,n?{4,-4},mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
id vert(+6,n?{4,-4},-6,mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
id vert(-6,n?{4,-4},+6,mu?,al?,be?,p?,q?,k?) =               -1/2*g*pm(4)*pm(4)/pm(3);
* ----------------------------------
* --- 3 4 6 --- W h phi+- ---
id vert(+3,-6,n?{4,-4},mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (q(mu)-k(mu));
id vert(-3,+6,n?{4,-4},mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (q(mu)-k(mu));
id vert(+3,n?{4,-4},-6,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (k(mu)-q(mu));
id vert(-3,n?{4,-4},+6,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (k(mu)-q(mu));
id vert(n?{4,-4},+3,-6,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (k(al)-p(al));
id vert(n?{4,-4},-3,+6,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (k(al)-p(al));
id vert(n?{4,-4},+6,-3,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (q(be)-p(be));
id vert(n?{4,-4},-6,+3,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (q(be)-p(be));
id vert(+6,-3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (p(al)-k(al));
id vert(-6,+3,n?{4,-4},mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (p(al)-k(al));
id vert(+6,n?{4,-4},-3,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (p(be)-q(be));
id vert(-6,n?{4,-4},+3,mu?,al?,be?,p?,q?,k?) =                i_*g/2*    (p(be)-q(be));
* --- 2 4 5 --- Z h phi0 ---
id vert(n?{2,-2},n2?{5,-5},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(q(mu)-k(mu));
id vert(n?{2,-2},n2?{4,-4},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(k(mu)-q(mu));
id vert(n?{4,-4},n2?{2,-2},n3?{5,-5},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(k(al)-p(al));
id vert(n?{4,-4},n2?{5,-5},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(q(be)-p(be));
id vert(n?{5,-5},n2?{2,-2},n3?{4,-4},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(p(al)-k(al));
id vert(n?{5,-5},n2?{4,-4},n3?{2,-2},mu?,al?,be?,p?,q?,k?) = +i_*g/2/ctw*(p(be)-q(be));
* --- 3 5 6 --- W phi0 phi+- ---
id vert(+3,-6,n?{5,-5},mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (q(mu)-k(mu));
id vert(-3,+6,n?{5,-5},mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (k(mu)-q(mu));
id vert(+3,n?{5,-5},-6,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (k(mu)-q(mu));
id vert(-3,n?{5,-5},+6,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (q(mu)-k(mu));
id vert(n?{5,-5},+3,-6,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (k(al)-p(al));
id vert(n?{5,-5},-3,+6,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (p(al)-k(al));
id vert(n?{5,-5},+6,-3,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (p(be)-q(be));
id vert(n?{5,-5},-6,+3,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (q(be)-p(be));
id vert(+6,-3,n?{5,-5},mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (k(al)-p(al));
id vert(-6,+3,n?{5,-5},mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (p(al)-k(al));
id vert(+6,n?{5,-5},-3,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (q(be)-p(be));
id vert(-6,n?{5,-5},+3,mu?,al?,be?,p?,q?,k?) =                 1/2*g*    (p(be)-q(be));
* --- 1 6 6 --- gamma phi+- phi+- ---
id vert(n?{1,-1},-6,+6,mu?,al?,be?,p?,q?,k?) =        +g*stw*            (k(mu)-q(mu));
id vert(n?{1,-1},+6,-6,mu?,al?,be?,p?,q?,k?) =        +g*stw*            (q(mu)-k(mu));
id vert(+6,-6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =        +g*stw*            (p(be)-q(be));
id vert(-6,+6,n?{1,-1},mu?,al?,be?,p?,q?,k?) =        +g*stw*            (q(be)-p(be));
id vert(+6,n?{1,-1},-6,mu?,al?,be?,p?,q?,k?) =        +g*stw*            (p(al)-k(al));
id vert(-6,n?{1,-1},+6,mu?,al?,be?,p?,q?,k?) =        +g*stw*            (k(al)-p(al));
* --- 2 6 6 --- Z phi+ phi- ---                       
id vert(n?{2,-2},-6,+6,mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(k(mu)-q(mu));
id vert(n?{2,-2},+6,-6,mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(q(mu)-k(mu));
id vert(+6,-6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(p(be)-q(be));
id vert(-6,+6,n?{2,-2},mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(q(be)-p(be));
id vert(+6,n?{2,-2},-6,mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(p(al)-k(al));
id vert(-6,n?{2,-2},+6,mu?,al?,be?,p?,q?,k?) =        +g*(2*ctw-1/ctw)/2*(k(al)-p(al));
* ----------------------------------
* --- 6 7 10 --- phi+- Xp Ya ---               
id vert(+6,-7,+10,mu?,al?,be?,p?,q?,k?)=       -i_*g*stw*rxi(3)*pm(3);
id vert(+6,+10,-7,mu?,al?,be?,p?,q?,k?)=       -i_*g*stw*rxi(3)*pm(3);
id vert(-6,-8,+10,mu?,al?,be?,p?,q?,k?)=        i_*g*stw*rxi(3)*pm(3);
id vert(-6,+10,-8,mu?,al?,be?,p?,q?,k?)=        i_*g*ctw*rxi(3)*pm(3);
* --- 6 7 9 --- phi+- Xp Yz ---
id vert(-6,-9,+7,mu?,al?,be?,p?,q?,k?) =       -i_/2*g/ctw*rxi(2)*pm(3);
id vert(-6,+7,-9,mu?,al?,be?,p?,q?,k?) =       -i_/2*g/ctw*rxi(2)*pm(3);
id vert(+6,-7,+9,mu?,al?,be?,p?,q?,k?) =       -i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
id vert(+6,+9,-7,mu?,al?,be?,p?,q?,k?) =       -i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
* --- 6 8 9 --- phi+- Xm Yz ---
id vert(+6,-9,+8,mu?,al?,be?,p?,q?,k?) =        i_/2*g/ctw*rxi(2)*pm(3);
id vert(+6,+8,-9,mu?,al?,be?,p?,q?,k?) =        i_/2*g/ctw*rxi(2)*pm(3);
id vert(-6,-8,+9,mu?,al?,be?,p?,q?,k?) =        i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);
id vert(-6,+9,-8,mu?,al?,be?,p?,q?,k?) =        i_*g*(2*ctw-1/ctw)/2*rxi(3)*pm(3);

* --- 3 7 9 --- W Xp Yz ---
id vert(+3,-7,+9,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(3)*q(mu);
id vert(+3,+9,-7,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(3)*k(mu);
id vert(-3,-9,+7,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(2)*q(mu);
id vert(-3,+7,-9,mu?,al?,be?,p?,q?,k?)  =          +g*ctw/rxi(2)*k(mu);
* --- 3 7 10 --- W Xp Ya ---
id vert(+3,-7,+10,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(3)*q(mu);
id vert(+3,+10,-7,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(3)*k(mu);
id vert(-3,-10,+7,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(1)*q(mu);
id vert(-3,+7,-10,mu?,al?,be?,p?,q?,k?) =          +g*stw/rxi(1)*k(mu);
* --- 3 8 9 --- W Xm Yz ---
id vert(+3,-9,+8,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(2)*q(mu);
id vert(+3,+8,-9,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(2)*k(mu);
id vert(-3,-8,+9,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(3)*q(mu);
id vert(-3,+9,-8,mu?,al?,be?,p?,q?,k?)  =          -g*ctw/rxi(3)*k(mu);
* --- 3 8 10 --- W Xm Ya ---
id vert(+3,-10,+8,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(1)*q(mu);
id vert(+3,+8,-10,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(1)*k(mu);
id vert(-3,-8,+10,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(3)*q(mu);
id vert(-3,+10,-8,mu?,al?,be?,p?,q?,k?) =          -g*stw/rxi(3)*k(mu);

* --- 1 7 7 --- gamma Xp Xp ---
id vert(n?{1,-1},-7,+7,mu?,al?,be?,p?,q?,k?) =     -g*stw*q(mu)/rxi(3);
id vert(n?{1,-1},+7,-7,mu?,al?,be?,p?,q?,k?) =     -g*stw*k(mu)/rxi(3);
* --- 2 7 7 --- Z Xp Xp ---
id vert(n?{2,-2},-7,+7,mu?,al?,be?,p?,q?,k?) =     -g*ctw*q(mu)/rxi(3);
id vert(n?{2,-2},+7,-7,mu?,al?,be?,p?,q?,k?) =     -g*ctw*k(mu)/rxi(3);
* --- 4 7 7 --- h Xp Xp ---
id vert(n?{4,-4},+7,-7,mu?,al?,be?,p?,q?,k?) =     -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},-7,+7,mu?,al?,be?,p?,q?,k?) =     -1/2*g*pm(3)*rxi(3);
* --- 5 7 7 --- phi0 Xp Xp ---
id vert(n?{5,-5},+7,-7,mu?,al?,be?,p?,q?,k?) =     i_/2*g*pm(3)*rxi(3);
id vert(n?{5,-5},-7,+7,mu?,al?,be?,p?,q?,k?) =     i_/2*g*pm(3)*rxi(3);

* --- 1 8 8 --- gamma Xm Xm ---
id vert(n?{1,-1},-8,+8,mu?,al?,be?,p?,q?,k?) =     +g*stw*q(mu)/rxi(3);
id vert(n?{1,-1},+8,-8,mu?,al?,be?,p?,q?,k?) =     +g*stw*k(mu)/rxi(3);
* --- 2 8 8 --- Z Xm Xm ---
id vert(n?{2,-2},-8,+8,mu?,al?,be?,p?,q?,k?) =     +g*ctw*q(mu)/rxi(3);
id vert(n?{2,-2},+8,-8,mu?,al?,be?,p?,q?,k?) =     +g*ctw*k(mu)/rxi(3);
* --- 4 8 8 --- h Xm Xm ---
id vert(n?{4,-4},+8,-8,mu?,al?,be?,p?,q?,k?) =     -1/2*g*pm(3)*rxi(3);
id vert(n?{4,-4},-8,+8,mu?,al?,be?,p?,q?,k?) =     -1/2*g*pm(3)*rxi(3);
* --- 5 8 8 --- phi0 Xm Xm ---
id vert(n?{5,-5},+8,-8,mu?,al?,be?,p?,q?,k?) =     -i_/2*g*pm(3)*rxi(3);
id vert(n?{5,-5},-8,+8,mu?,al?,be?,p?,q?,k?) =     -i_/2*g*pm(3)*rxi(3);

* --- 4 9 9 --- h Yz Yz ---
id vert(n?{4,-4},+9,-9,mu?,al?,be?,p?,q?,k?) =     -1/2*g/ctw/ctw*pm(3)*rxi(2);
id vert(n?{4,-4},-9,+9,mu?,al?,be?,p?,q?,k?) =     -1/2*g/ctw/ctw*pm(3)*rxi(2);


*           ---------boson_boson_four_line_vertices-----------
*
* --- 1 1 3 3 --- gamma gamma Z Z --- (12 Rules)
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
* --- 2 2 3 3 --- Z Z W W --- (12 Rules)
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
* --- 1 2 3 3 --- gamma Z W W --- (24 Rules)
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
* --- 3 3 3 3 --- W W W W ---
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

* --- 1 1 6 6 --- gamma gamma phi+ phi- ---
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
* --- 2 2 6 6 --- Z Z phi+ phi- ---
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
* --- 1 2 6 6 --- gamma Z phi+ phi- ---
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

* --- 1 3 4 6 --- gamma W h phi+- ---
id vert(n1?{4,-4},n2?{1,-1},+6,-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},n2?{1,-1},-6,+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},n2?{1,-1},+3,-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(n1?{4,-4},n2?{1,-1},-3,+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
*
id vert(n1?{4,-4},+6,n2?{1,-1},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},-6,n2?{1,-1},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},+3,n2?{1,-1},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(n1?{4,-4},-3,n2?{1,-1},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
*
id vert(n1?{4,-4},+6,-3,n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},-6,+3,n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(n1?{4,-4},+3,-6,n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(n1?{4,-4},-3,+6,n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
*
id vert(+6,-3,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
id vert(-6,+3,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(+3,-6,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,be);
id vert(-3,+6,n1?{4,-4},n2?{1,-1},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,be);
*
id vert(+6,-3,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,al);
id vert(-6,+3,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,al);
id vert(+3,-6,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(-3,+6,n1?{1,-1},n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
*
id vert(n1?{1,-1},+6,-3,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},-6,+3,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},+3,-6,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,nu);
id vert(n1?{1,-1},-3,+6,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,nu);
*
id vert(-3,n1?{4,-4},n2?{1,-1},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,al);
id vert(+3,n1?{4,-4},n2?{1,-1},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,al);
id vert(-6,n1?{4,-4},n2?{1,-1},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(al,be);
id vert(+6,n1?{4,-4},n2?{1,-1},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(al,be);
*
id vert(-3,n2?{1,-1},n1?{4,-4},+6,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(mu,nu);
id vert(+3,n2?{1,-1},n1?{4,-4},-6,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(mu,nu);
id vert(-6,n2?{1,-1},n1?{4,-4},+3,mu?,nu?,al?,be?)=  i_/2*g^2*stw*d_(nu,be);
id vert(+6,n2?{1,-1},n1?{4,-4},-3,mu?,nu?,al?,be?)= -i_/2*g^2*stw*d_(nu,be);
* --- 1 3 5 6 --- gamma Z phi0 phi+- ---
id vert(+6,-3,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(-6,+3,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(+3,-6,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(-3,+6,n1?{5,-5},n2?{1,-1},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
*
id vert(+6,-3,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(-6,+3,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(+3,-6,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
id vert(-3,+6,n1?{1,-1},n2?{5,-5},mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
*
id vert(n1?{5,-5},n2?{1,-1},+6,-3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(n1?{5,-5},n2?{1,-1},-6,+3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,be);
id vert(n1?{5,-5},n2?{1,-1},+3,-6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
id vert(n1?{5,-5},n2?{1,-1},-3,+6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(nu,al);
*
id vert(n1?{1,-1},n2?{5,-5},+6,-3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(n1?{1,-1},n2?{5,-5},-6,+3,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,be);
id vert(n1?{1,-1},n2?{5,-5},+3,-6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
id vert(n1?{1,-1},n2?{5,-5},-3,+6,mu?,nu?,al?,be?)=  1/2*g^2*stw*d_(mu,al);
* --- 2 3 5 6 --- Z W phi0 phi+- ---
id vert(n1?{5,-5},n2?{2,-2},+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{5,-5},n2?{2,-2},-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,be);
id vert(n1?{5,-5},n2?{2,-2},+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,al);
id vert(n1?{5,-5},n2?{2,-2},-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(nu,al);
*
id vert(n1?{2,-2},n2?{5,-5},+6,-3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},-6,+3,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},+3,-6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},n2?{5,-5},-3,+6,mu?,nu?,al?,be?)= -1/2*g^2*stw^2/ctw*d_(mu,al);
*-2634 dobavlyat'....
id vert(n1?{2,-2},+6,-3,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},-6,+3,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,al);
id vert(n1?{2,-2},+3,-6,n2?{4,-4},mu?,nu?,al?,be?)= -i_/2*g^2*stw^2/ctw*d_(mu,nu);
id vert(n1?{2,-2},-3,+6,n2?{4,-4},mu?,nu?,al?,be?)=  i_/2*g^2*stw^2/ctw*d_(mu,nu);
* --- 2 3 4 6 --- Z W h phi+- ---
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

* --- 2 2 4 4 --- Z Z h h ---
id vert(n1?{2,-2},n2?{2,-2},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,nu);
id vert(n1?{2,-2},n2?{4,-4},n3?{4,-4},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,be);
id vert(n1?{2,-2},n2?{4,-4},n3?{2,-2},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,al);
id vert(n1?{4,-4},n2?{2,-2},n3?{2,-2},n4?{4,-4},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,al);
id vert(n1?{4,-4},n2?{2,-2},n3?{4,-4},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,be);
id vert(n1?{4,-4},n2?{4,-4},n3?{2,-2},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw/ctw*d_(al,be);
* --- 2 2 5 5 --- Z Z phi0 phi0 ---
id vert(n1?{2,-2},n2?{2,-2},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,nu);
id vert(n1?{2,-2},n2?{5,-5},n3?{5,-5},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,be);
id vert(n1?{2,-2},n2?{5,-5},n3?{2,-2},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(mu,al);
id vert(n1?{5,-5},n2?{2,-2},n3?{2,-2},n4?{5,-5},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,al);
id vert(n1?{5,-5},n2?{2,-2},n3?{5,-5},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(nu,be);
id vert(n1?{5,-5},n2?{5,-5},n3?{2,-2},n4?{2,-2},mu?,nu?,al?,be?)= -1/2*g^2/ctw^2*d_(al,be);
* --- 3 3 4 4 --- W W h h ---
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
* --- 3 3 5 5 --- W W phi0 phi0 ---
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
* --- 3 3 6 6 --- W W phi+ phi- ---
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

* --- 4 4 4 4 --- h h h h ---
id vert(n1?{4,-4},n2?{4,-4},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)=-3/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
* --- 5 5 5 5 --- phi0 phi0 phi0 phi0 ---
id vert(n1?{5,-5},n2?{5,-5},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)=-3/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
* --- 4 4 5 5 --- phi0 phi0 phi+ phi- ---
id vert(n1?{4,-4},n2?{4,-4},n3?{5,-5},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},n2?{5,-5},n3?{5,-5},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{4,-4},n2?{5,-5},n3?{4,-4},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{4,-4},n3?{5,-5},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{4,-4},n3?{4,-4},n4?{5,-5},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(n1?{5,-5},n2?{5,-5},n3?{4,-4},n4?{4,-4},mu?,nu?,al?,be?)=-1/4*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
* --- 4 4 6 6 --- h h phi+ phi- ---
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
* --- 5 5 6 6 --- phi0 phi0 phi+ phi- ---
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
* --- 6 6 6 6 --- phi+ phi- phi + phi- ---
id vert(+6,-6,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,-6,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,+6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,+6,-6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(-6,-6,+6,+6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;
id vert(+6,+6,-6,-6,mu?,nu?,al?,be?)= -1/2*g^2*pm(4)*pm(4)*pm(3)^-1*pm(3)^-1;

*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             11= nu_e      en     12= electron  el                  *
*             13= up        up     14= down      dn                  *
*             15= nu_mu     mn     16= muon      mu                  *
*             17= charm     ch     18= strange   st                  *
*             19= nu_tau    tn     20= tauon     ta                  *
*             21= top       tp     22= bottom    bt                  *
*====================================================================*
*
*-             ----------boson_fermion_vertices-----------
*
* first generation
* ----------------
* QCD rules for boson-quark vertices
id vert(n?{1,-1},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{2,-2},+13,-13,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-13,+13,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+14,-14,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-14,+14,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(-3,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(13)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(13)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(14)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(14)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{5,-5},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(13)-vma(13))/4/ctw*pm(13)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(14)-vma(14))/4/ctw*pm(14)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(14)*pm(3)^-1*gd6(ii)-pm(13)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(13)*pm(3)^-1*gd6(ii)-pm(14)*pm(3)^-1*gd7(ii))*d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+11,-11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-11,+11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+12,-12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-12,+12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+13,-13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-13,+13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+14,-14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-14,+14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(2);
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
id vert(n?{4,-4},+11,-11,mu?,ii?) = -g/4*pm(11)*pm(3)^-1*(2);
id vert(n?{4,-4},-11,+11,mu?,ii?) = -g/4*pm(11)*pm(3)^-1*(2);
id vert(n?{4,-4},+12,-12,mu?,ii?) = -g/4*pm(12)*pm(3)^-1*(2);
id vert(n?{4,-4},-12,+12,mu?,ii?) = -g/4*pm(12)*pm(3)^-1*(2);
id vert(n?{4,-4},+13,-13,mu?,ii?) = -g/4*pm(13)*pm(3)^-1*(2);
id vert(n?{4,-4},-13,+13,mu?,ii?) = -g/4*pm(13)*pm(3)^-1*(2);
id vert(n?{4,-4},+14,-14,mu?,ii?) = -g/4*pm(14)*pm(3)^-1*(2);
id vert(n?{4,-4},-14,+14,mu?,ii?) = -g/4*pm(14)*pm(3)^-1*(2);
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
id vert(n?{1,-1},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{2,-2},+17,-17,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-17,+17,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+18,-18,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-18,+18,mu?,cl1?scl,cl2?scl,ii?) =   i_*g/4/ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(-3,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(17)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(17)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(18)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(18)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{5,-5},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(17)-vma(17))/4/ctw*pm(17)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(18)-vma(18))/4/ctw*pm(18)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(18)*pm(3)^-1*gd6(ii)-pm(17)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(17)*pm(3)^-1*gd6(ii)-pm(18)*pm(3)^-1*gd7(ii))*d(cl1,cl2);

* SM vertices
id vert(n?{1,-1},+15,-15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-15,+15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+16,-16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-16,+16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+17,-17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-17,+17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+18,-18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-18,+18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(2);
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
id vert(n?{4,-4},+15,-15,mu?,ii?) = -g/4*pm(15)*pm(3)^-1*(2);
id vert(n?{4,-4},-15,+15,mu?,ii?) = -g/4*pm(15)*pm(3)^-1*(2);
id vert(n?{4,-4},+16,-16,mu?,ii?) = -g/4*pm(16)*pm(3)^-1*(2);
id vert(n?{4,-4},-16,+16,mu?,ii?) = -g/4*pm(16)*pm(3)^-1*(2);
id vert(n?{4,-4},+17,-17,mu?,ii?) = -g/4*pm(17)*pm(3)^-1*(2);
id vert(n?{4,-4},-17,+17,mu?,ii?) = -g/4*pm(17)*pm(3)^-1*(2);
id vert(n?{4,-4},+18,-18,mu?,ii?) = -g/4*pm(18)*pm(3)^-1*(2);
id vert(n?{4,-4},-18,+18,mu?,ii?) = -g/4*pm(18)*pm(3)^-1*(2);
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
id vert(n?{1,-1},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(n1)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{1,-1},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(n2)*stw*gd(ii,mu)*(2)*d(cl1,cl2);
id vert(n?{2,-2},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4/ctw*gd(ii,mu)*(vpa(n1)*gd6(ii)+vma(n1)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4/ctw*gd(ii,mu)*(vpa(n2)*gd6(ii)+vma(n2)*gd7(ii))*d(cl1,cl2);
id vert(-3,+21,-22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(-3,-22,+21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,+22,-21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(+3,-21,+22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*gd(ii,mu)*gd6(ii)*d(cl1,cl2);
id vert(n?{4,-4},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(n1)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{4,-4},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = -g/4*pm(n2)*pm(3)^-1*(2)*d(cl1,cl2);
id vert(n?{5,-5},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(n1)-vma(n1))/4/ctw*pm(n1)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(n?{5,-5},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g*(vpa(n2)-vma(n2))/4/ctw*pm(n2)*pm(2)^-1*(gd6(ii)-gd7(ii))*d(cl1,cl2);
id vert(-6,+21,-22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(-6,-22,+21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(22)*pm(3)^-1*gd6(ii)-pm(21)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,+22,-21,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
id vert(+6,-21,+22,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2/sr2*(pm(21)*pm(3)^-1*gd6(ii)-pm(22)*pm(3)^-1*gd7(ii))*d(cl1,cl2);
* SM vertices
id vert(n?{1,-1},+19,-19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-19,+19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+20,-20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-20,+20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+21,-21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-21,+21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},+22,-22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(2);
id vert(n?{1,-1},-22,+22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(2);
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
id vert(n?{4,-4},+19,-19,mu?,ii?) = -g/4*pm(19)*pm(3)^-1*(2);
id vert(n?{4,-4},-19,+19,mu?,ii?) = -g/4*pm(19)*pm(3)^-1*(2);
id vert(n?{4,-4},+20,-20,mu?,ii?) = -g/4*pm(20)*pm(3)^-1*(2);
id vert(n?{4,-4},-20,+20,mu?,ii?) = -g/4*pm(20)*pm(3)^-1*(2);
id vert(n?{4,-4},+21,-21,mu?,ii?) = -g/4*pm(21)*pm(3)^-1*(2);
id vert(n?{4,-4},-21,+21,mu?,ii?) = -g/4*pm(21)*pm(3)^-1*(2);
id vert(n?{4,-4},+22,-22,mu?,ii?) = -g/4*pm(22)*pm(3)^-1*(2);
id vert(n?{4,-4},-22,+22,mu?,ii?) = -g/4*pm(22)*pm(3)^-1*(2);
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
*
id vert(?a)=0;
id sr2^-2 = 1/2;

* ###################################################
* #            P R O P A G A T O R S                #
* ###################################################
*
*   
* ========== SM Bosons ==========
* ===============================
*
*id pr(n?{1,-1},mu?,nu?,p?)  = den(1,0,p)*(d_(mu,nu) + (rxi(1)^2-1)*p(mu)*p(nu)*den(1,0,p));
*id pr(n?{2,-2},mu?,nu?,p?)  = den(1,pm(2),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(2)^-1*pm(2)^-1)
*                              - p(mu)*p(nu)*pm(2)^-1*pm(2)^-1*den(rxi(2),pm(2),p);
*id pr(n?{3,-3},mu?,nu?,p?)  = den(1,pm(3),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(3)^-1*pm(3)^-1)
*                              - p(mu)*p(nu)*pm(3)^-1*pm(3)^-1*den(rxi(3),pm(3),p);

id pr(n?{1,-1},mu?,nu?,p?)  = den(1,0,p)*(d_(mu,nu));
id pr(n?{2,-2},mu?,nu?,p?)  = den(1,pm(2),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(2)^-1*pm(2)^-1);
id pr(n?{3,-3},mu?,nu?,p?)  = den(1,pm(3),p)*(d_(mu,nu) + p(mu)*p(nu)*pm(3)^-1*pm(3)^-1);
id pr(n?{4,-4},mu?,nu?,p?)  = den(1,pm(4),p);

* == Fermions - Quarks/leptons ==
* ===============================
*
id pr(n?sQI,?a,cl1?scl,cl2?scl,i?) = d(cl1,cl2)*den(1,pm(abs_(n)),?a)*pr(sf,n,?a,i);
id pr(n?sFI,p?,i?) = den(1,pm(abs_(n)),p)*pr(sf,n,p,i);

splitarg pr 3;
#if `PeskinNaumov'
 repeat id pr(sf,n?sFI,p?,?a,i?) = i_*gd(i,p)+pr(sf,n,?a,i); 
 id pr(sf,n?sFI,i?) = i_*pm(abs_(n));
#else
 repeat id pr(sf,n?sFI,p?,?a,i?) = -i_*gd(i,p)+pr(sf,n,?a,i); 
 id pr(sf,n?sFI,i?) = pm(abs_(n));
#endif
id gd(i?,-p?) = -gd(i,p);


 repeat;
   id d(?a,cl1?sAcl,?b)*d(?c,cl1?sAcl,?d) = d(?a,?c,?d,?b);
   id d(?a,cl1?sAgi,?b)*d(?c,cl1?sAgi,?d) = d(?a,?c,?d,?b);
 endrepeat;

* I=0: QED,QCD
* I=1: EW
#if {`I'} = 0
    #call a2b(gd6,gd5)
    #call a2b(gd7,gd5)
    #call a2b(g,e)
#endif

#endprocedure
