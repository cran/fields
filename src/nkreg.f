
c*** computes a  regression estimate for the multivariate data xk
c**** with bandwidth h. The desnity function is evaluated at
c**** the points x
c**** kernel used is a  mulitvariate normal
c**** 
c**** 
      subroutine nkreg( h,n,m,xk,y,np,x,f)
      REAL*8 xk(n,m),y(n),H,f(n),x(np,m)
      real*8 tempd, tempn,tempw, dist, rh2
      integer n,m,np,j,i,jj

      rh2= 1.0/(h*h)
      DO 2 J=1,NP              

      tempd= 0.0
      tempn= 0.0

         DO 5 I= 1,n                                                  

          dist=0.0
c**** loop over components of the vectors 
          do 6 jj=1,m
              dist= dist +  ( x(j,jj)-xk(i,jj))**2
  6       continue 
c**** accumulate contribution from each data point
          tempw= exp( -.5*(dist)*rh2)
          tempn= tempn + tempw* y(i) 
          tempd=tempd + tempw
          
  5     continue

c**** normalize estimate
      f(j)= tempn/tempd


c**** f(j) = estimate of  regression function at the point [x(j,1),...,x(j,m)]

  2   continue                                   
      return    
      END    
