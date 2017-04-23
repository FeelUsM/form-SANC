#procedure Stop(message)
***********************************
* Terminates the program after    *
* execution of a current module   *
* has finished.                   *
*                                 *
* message -- printed before exit  * 
***********************************

.sort
 #if (`$error' == 1)        
      #write " ## `message'";
     .end
 #endif     
#endprocedure
*------------
