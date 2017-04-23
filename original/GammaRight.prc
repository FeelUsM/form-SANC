#procedure GammaRight(dummy)
*---------------------------
 repeat;
   id gd5(i?)*gd5(i?) = 1;
   id gd6(i?)*gd6(i?) = 2*gd6(i);
   id gd7(i?)*gd7(i?) = 2*gd7(i);
   id gd6(i?)*gd7(i?) = 0;
   id gd7(i?)*gd6(i?) = 0;
   id gd5(i?)*gd6(i?) = +gd6(i);
   id gd5(i?)*gd7(i?) = -gd7(i);
   id gd7(i?)*gd5(i?) = -gd7(i);
   id gd6(i?)*gd5(i?) = +gd6(i);
   id gd5(i?)*gd(i?,mu1?) = -gd(i,mu1)*gd5(i);
   id gd6(i?)*gd(i?,mu1?) = +gd(i,mu1)*gd7(i);
   id gd7(i?)*gd(i?,mu1?) = +gd(i,mu1)*gd6(i);
 endrepeat;

#endprocedure
*------------ 
