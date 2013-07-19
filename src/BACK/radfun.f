
       subroutine radfun(n,d2, par)
       real*8 d2(n), par(2), dtemp
       integer n
       if( int(par(2)).eq.0) then
         
         do 5 k =1,n
           dtemp= d2(k)
           if( dtemp.lt.1e-20) dtemp =1e-20 
         d2(k)= (dtemp)**( par(1))
   5     continue
        else 
         do 6 k=1,n
          dtemp= d2(k)
          if( dtemp.gt.1e-20)  then
c note: dtemp is squared distance
c divide by 2 to have log evaluated on distance
c as opposed to squared distance. 

           d2(k)= log(dtemp)*(dtemp)**( par(1))/2
          else
           d2(k)=0.0
          endif
   6   continue
       endif
        return
        end

