#if `FeynmanH'
#redefine FeynmanH "0"

* === Дираковские матрицы ===
* Nfun	axgd,gd4,gdI,cR,pr;				* Declar.h:283          * gd4 - Dirac gamma_{4}...
Nfun	gd,gd5,gd6,gd7;					* Declar.h:260
Set sgd:gd,gd5,gd6,gd7;                   * set of Dirac Gamma matrices... 
* Set sbg:g_,g5_,g6_,g7_;                   * set of Built-in Dirac Gamma matrices... 

* === константы электрослабого взаимодействия ===
sym	e,g,stw,ctw,GFermi,sr2; *,gs,CFScheme,qW;
* sr2 - какой-то коэф-т в W-бозонных вершинах

* === электрические заряды ===
sym	qel,qup,qdn;
sym	qmo,qch,qst; 
sym	qta,qtp,qbt;

* === коэффициенты левости и правости ===
sym	vpaen,vmaen,vpael,vmael,vpaup,vmaup,vpadn,vmadn;
sym	vpamn,vmamn,vpamo,vmamo,vpach,vmach,vpast,vmast;
sym	vpatn,vmatn,vpata,vmata,vpatp,vmatp,vpabt,vmabt;
* = еще что-то с этим связанное =
sym	ven,vel,vup,vdn, aen,ael,aup,adn, i3en,i3el,i3up,i3dn;
sym	vmn,vmo,vch,vst, amn,amo,ach,ast, i3mn,i3mo,i3ch,i3st;
sym	vtn,vta,vtp,vbt, atn,ata,atp,abt, i3tn,i3ta,i3tp,i3bt;

* === delta-function for QCD ===
cfun	d; * d(color1,color2)

* === denumerator ===
cfun	den; * den(1,m,p) = 1/(p.p-m^2)

#endif
