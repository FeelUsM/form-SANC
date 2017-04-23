#procedure SetFlags()
******************/
* Empty realisation of SetFlags   *
* procedure. One should overload  *
* this procedure inside program.  *
* Some procedures uses this one,  *
* which speeds up the calculation.*
******************/
 #ifdef `debug'
    #write " ## >>>>>> SetFlags() - EMPTY PROCEDURE: ";
    #write "           redefine inside your program without .sort and .store !";
 #endif
#endprocedure
*------------
