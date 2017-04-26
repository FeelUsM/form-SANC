#procedure Feynman(I)
*----------------------------
* I=0: QED,QCD
* I=1: EW

#define gd6pgd7eq2 "0"
* устраняется при помощи #call mya2b(gd7,gd6) или еще проще Feynman(0)

.sort :FeynmanStart;
#include Misc.h
vector p,q,k; * в трехбозонных вершинах
Sym n1,n2;
index i;
Cfun	sf;	* for internal use
#include myDeclare.h
#include Feynman.h;
#call myGlobal()



* ###########################################################################################
* #                                                                                         #
* #                                FEYNMAN RULES FOR SM                                     #
* #                                                                                         #
* ###########################################################################################

*====================================================================*
*                     ------- List of fields ------------            *
*                                                                    *
*             1 = gamma            2 = z           +/-3 = w^{+/-}    *
*====================================================================*
*             -------------boson_boson_vertices-----------
* === gamma,w+,w- ===
id vert(n?{1,-1},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{1,-1},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*stw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
* === z,w+,w- ===
id vert(n?{2,-2},+3,-3,mu?,al?,be?,p?,q?,k?) =    +g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
id vert(n?{2,-2},-3,+3,mu?,al?,be?,p?,q?,k?) =    -g*ctw*(d_(mu,al)*(p(be)-q(be))
                                                         +d_(al,be)*(q(mu)-k(mu))
                                                         +d_(mu,be)*(k(al)-p(al)));
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
#if `gd6pgd7eq2'
id vert(n?{1,-1},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*2*d(cl1,cl2);
#else
id vert(n?{1,-1},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
#endif
id vert(n?{2,-2},+13,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-13,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+14,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-14,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+13,-14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-14,+13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+14,-13,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-13,+14,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);

* SM vertices
#if `gd6pgd7eq2'
id vert(n?{1,-1},+11,-11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-11,+11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+12,-12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-12,+12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+13,-13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-13,+13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+14,-14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-14,+14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*2;
#else
id vert(n?{1,-1},+11,-11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-11,+11,mu?,ii?) = i_/2*g*c(11)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+12,-12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-12,+12,mu?,ii?) = i_/2*g*c(12)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+13,-13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-13,+13,mu?,ii?) = i_/2*g*c(13)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+14,-14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-14,+14,mu?,ii?) = i_/2*g*c(14)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
#endif
id vert(n?{2,-2},+11,-11,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},-11,+11,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(11)*gd6(ii)+vma(11)*gd7(ii));
id vert(n?{2,-2},+12,-12,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},-12,+12,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(12)*gd6(ii)+vma(12)*gd7(ii));
id vert(n?{2,-2},+13,-13,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},-13,+13,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(13)*gd6(ii)+vma(13)*gd7(ii));
id vert(n?{2,-2},+14,-14,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
id vert(n?{2,-2},-14,+14,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(14)*gd6(ii)+vma(14)*gd7(ii));
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
* second generation
* -----------------
* QCD rules for boson-quark vertices
#if `gd6pgd7eq2'
id vert(n?{1,-1},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*2*d(cl1,cl2);
#else
id vert(n?{1,-1},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
#endif
id vert(n?{2,-2},+17,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-17,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},+18,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},-18,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+17,-18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-18,+17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+18,-17,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-17,+18,mu?,cl1?scl,cl2?scl,ii?) = i_/2*g      /sr2*gd(ii,mu)*         gd6(ii)                 *d(cl1,cl2);

* SM vertices
#if `gd6pgd7eq2'
id vert(n?{1,-1},+15,-15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-15,+15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+16,-16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-16,+16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+17,-17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-17,+17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+18,-18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-18,+18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*2;
#else
id vert(n?{1,-1},+15,-15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-15,+15,mu?,ii?) = i_/2*g*c(15)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+16,-16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-16,+16,mu?,ii?) = i_/2*g*c(16)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+17,-17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-17,+17,mu?,ii?) = i_/2*g*c(17)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+18,-18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-18,+18,mu?,ii?) = i_/2*g*c(18)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
#endif
id vert(n?{2,-2},+15,-15,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},-15,+15,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(15)*gd6(ii)+vma(15)*gd7(ii));
id vert(n?{2,-2},+16,-16,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},-16,+16,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(16)*gd6(ii)+vma(16)*gd7(ii));
id vert(n?{2,-2},+17,-17,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},-17,+17,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(17)*gd6(ii)+vma(17)*gd7(ii));
id vert(n?{2,-2},+18,-18,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
id vert(n?{2,-2},-18,+18,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(18)*gd6(ii)+vma(18)*gd7(ii));
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
* third generation
* ----------------
* QCD rules for boson-quark vertices
#if `gd6pgd7eq2'
id vert(n?{1,-1},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n1)*stw*gd(ii,mu)*2*d(cl1,cl2);
id vert(n?{1,-1},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n2)*stw*gd(ii,mu)*2*d(cl1,cl2);
#else
id vert(n?{1,-1},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n1)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
id vert(n?{1,-1},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2*c(n2)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii))*d(cl1,cl2);
#endif
id vert(n?{2,-2},n1?sIQI,n2?sOQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4      /ctw*gd(ii,mu)*(vpa(n1)*gd6(ii)+vma(n1)*gd7(ii))*d(cl1,cl2);
id vert(n?{2,-2},n1?sOQI,n2?sIQI,mu?,cl1?scl,cl2?scl,ii?) = i_*g/4      /ctw*gd(ii,mu)*(vpa(n2)*gd6(ii)+vma(n2)*gd7(ii))*d(cl1,cl2);
id vert(-3      ,+21    ,-22    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(-3      ,-22    ,+21    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,+22    ,-21    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
id vert(+3      ,-21    ,+22    ,mu?,cl1?scl,cl2?scl,ii?) = i_*g/2      /sr2*gd(ii,mu)         *gd6(ii)                 *d(cl1,cl2);
* SM vertices
#if `gd6pgd7eq2'
id vert(n?{1,-1},+19,-19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-19,+19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+20,-20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-20,+20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+21,-21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-21,+21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},+22,-22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*2;
id vert(n?{1,-1},-22,+22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*2;
#else
id vert(n?{1,-1},+19,-19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-19,+19,mu?,ii?) = i_/2*g*c(19)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+20,-20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-20,+20,mu?,ii?) = i_/2*g*c(20)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+21,-21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-21,+21,mu?,ii?) = i_/2*g*c(21)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},+22,-22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
id vert(n?{1,-1},-22,+22,mu?,ii?) = i_/2*g*c(22)*stw*gd(ii,mu)*(        gd6(ii)+        gd7(ii));
#endif
id vert(n?{2,-2},+19,-19,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},-19,+19,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(19)*gd6(ii)+vma(19)*gd7(ii));
id vert(n?{2,-2},+20,-20,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},-20,+20,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(20)*gd6(ii)+vma(20)*gd7(ii));
id vert(n?{2,-2},+21,-21,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},-21,+21,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(21)*gd6(ii)+vma(21)*gd7(ii));
id vert(n?{2,-2},+22,-22,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
id vert(n?{2,-2},-22,+22,mu?,ii?) = i_/4*g      /ctw*gd(ii,mu)*(vpa(22)*gd6(ii)+vma(22)*gd7(ii));
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
    #call mya2b(gd6,gd5)
    #call mya2b(gd7,gd5)
    #call mya2b(g,e)
#endif

#endprocedure
