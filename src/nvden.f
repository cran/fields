
      subroutine nvden( h,n,m,xk,np,x,f)
      REAL*8 XK(n,m),H(n),f(n),x(np,m)
      real*8 tempd, hm,  dist, k0,con,rh2
      integer n,m,np,j,i,jj

c***** test for ordered first column in data
c***** send back negative density if in error
c     k0= KERNEL(0.0)
      k0=   .39894228
      con=   ( k0**m)

      DO 2 J=1,NP              

      tempd= 0.0

      DO 5 I= 1,n                                                  
      hm= h(i)**m
      rh2= 1.0/(h(i)*h(i))

          dist=0.0
c**** loop over components of the vectors 
          do 6 jj=1,m
              dist= dist +  ( x(j,jj)-xk(i,jj))**2
  6       continue 
c         dist=  dist/h
c**** accumulate contribution from each data point
          tempd=tempd +    con*exp( -.5*(dist)*rh2)/hm
c          write(*,*) jj,dist, tempd
          

  5     continue

c**** normalize estimate
c     f(j)= con*tempd/(n*hm)
      f(j)= tempd/n


c**** f(j) = estimate of probability density at the point [x(j,1),...,x(j,m)]

  2   continue                                   
      return    

      end
