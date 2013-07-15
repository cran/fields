C** evaluates the derivatives of the  radial basis functions 
c* and multiplies them by the matrix c
c**** and does the multplication  h= Kc
c****  K is n1Xn2
c****  h is n1Xn3
c***** c is n2xn3
 

       subroutine derrb( nd,x1,n1, x2,n2, par, c,n3,h,work)
       implicit double precision (a-h,o-z)
       integer nd,n1,n2,n3,ic, jc,j
       
       real*8 par(2),x1(n1,nd), x2(n2,nd), c(n2,n3), h(n1,n3),sum
       real*8 work( n1), ddot

c****** work aray must be dimensioned to size n1
c **** loop through columns of output matrix K
c*** outer most loop over columns of x1 and x2 should reduce swapping

       do 5 ir= 1, n1

c
 
c evaluate all basis functions at  x1(j,.)       
       do 10 j =1,n2
c
c  zero out sum accumulator
c
         sum=0.0
      do 15  ic=1,nd
c
c** accumulate squared differences
c 

            sum= sum+ (x1(ir,ic)- x2(j,ic))**2

 15             continue
        work(j)=sum
 10    continue

C**** evaluate squared distances  with basis functions. 

          call radfun( n2,work(1),par)
c
c***** now the dot products you have all been waiting for!
c
       do 30 jc=1,n3
          h(ir,jc)= ddot( n2, work(1), 1, c(1,jc),1)
30     continue
  
 5      continue

       return
       end
