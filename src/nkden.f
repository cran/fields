c*** computes a density estimate for the multivariate data xk
c**** with bandwidth h. The desnity function is evaluated at
c**** the points x
c**** kernel used is a  mulitvariate normal
c**** 
c**** 
      subroutine nkden( h,n,m,xk,np,x,f)
      REAL*8 XK(n,m),H,f(n),x(np,m)
      real*8 tempd, hm,  dist, k0,con,rh2
      integer n,m,np,j,i,jj

c***** test for ordered first column in data
c***** send back negative density if in error
      hm= h**m
      rh2= 1.0/(h*h)
      k0=   .39894228
      con=    k0**m

      DO 2 J=1,NP              

      tempd= 0.0

         DO 5 I= 1,n                                                  

          dist=0.0
c**** loop over components of the vectors 
          do 6 jj=1,m
              dist= dist +  ( x(j,jj)-xk(i,jj))**2
  6       continue 
c         dist=  dist/h
c**** accumulate contribution from each data point
          tempd=tempd +    exp( -.5*(dist)*rh2)
c          write(*,*) jj,dist, tempd
          

  5     continue

c**** normalize estimate
      f(j)= con*tempd/(n*hm)


c**** f(j) = estimate of probability density at the point [x(j,1),...,x(j,m)]

  2   continue                                   
      return    
      END    
