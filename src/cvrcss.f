
       subroutine cvrcss(
c    arguments passed through to rcss
     + n,x,y,wt,sy,diag,din,dout,  
c
c  arguments needed for cv search
     +  nstep,maxit,trmin, trmax, hmin,hmax,hopt,vopt,
     +     tropt,mxstep,tabout,
     +                 ierr)
      implicit double precision (a-h,o-z)
      real*8  tabout(mxstep,4)
      real*8  x(n),y(n),wt(n),sy(n) 
      real*8 diag(n), din(10), dout(10)
      real*8  hopt,vopt,tropt,trmin,trmax,hmin, hmax
      data tau,tausq/.6180339d0,.3819660d0/
**** coarse search in bandwidth h: this search feeds (hl,hr) to golden search
c compute range for coarse search   
      if (mxstep .lt. nstep) then
         ierr=10
         return
      endif


c********* begin block for finding hmin and hmax

      if( hmin.ge.hmax) then
C*** for any spline problem it is assumed that lambda=1 is
c*** computable and can be used as a start value
      hmax=0.0      
      do  1020 k=1,25
          
      call  cvrf(hmax,n,x,y,wt,sy,trace,diag,din,dout,cv,ierr)

           if( trace.le.trmin)  goto 1030
            hmax= hmax + 1 
 1020 continue
       ierr=31
 1030 continue      

      hmin=0.0      
      do  1000 k=1,25
          call  cvrf(hmin,n,x,y,wt,sy,trace,diag,din,dout,cv,ierr)
           if( trace.ge.trmax)  goto 1010
            hmin= hmin - 1
 
 1000 continue
       ierr=21
 1010 continue      
      endif

c********* end block for finding hmin and hmax

      hstep=(hmax-hmin)/(nstep-1)
      do 3 j=0,nstep-1
         h=hmax - j*hstep
         call cvrf(h,n,x,y,wt,sy,trace,diag,din,dout,cvh,ierr)  
         tabout(j+1,1)=h
         tabout(j+1,2)=dout(3)
         tabout(j+1,3)=cvh
         tabout(j+1,4)=dout(1)
         if ((cvh .lt. cvmin) .or. (j .eq. 0)) then
            hopt=h
            best=h
            cvmin=cvh
            trbest=trace
         endif
 5000   format(5e12.4)
c
  3   continue
c
c**** fast return if crude search minimum cv is hmin or hmax
           hopt= best
           vopt=cvmin
           tropt=trbest
      if( (best.le.hmin) .or. (best.ge.hmax) ) then
           ierr=-1
           return
      endif


c**** start values for golden search
      hl=best-hstep
      hr=best+hstep

*** Golden section search for min cv on interval (hl,hr)--maxit iterations
***   cv(h) must be quasiconvex on initial interval (hl,hr).
***   On return, h=abscissa of minimum, v=cv(h)=actual minimum achieved.
***   Interval of uncertainty is (hl,hr), hl <= hlm <= hrm <= hr.

      call 
     * cvrf(hl,n,x,y,wt,sy,trace,diag,din,dout,cvhl,ierr)  

      call 
     * cvrf(hr,n,x,y,wt,sy,trace,diag,din,dout,cvhr,ierr)  

      hlm=hl*tau+hr*tausq
      hrm=hl+hr-hlm
      call 
     * cvrf(hlm,n,x,y,wt,sy,trchlm,diag,din,dout,cvhlm,ierr)  

      call 
     * cvrf(hrm,n,x,y,wt,sy,trchrm,diag,din,dout,cvhrm,ierr)  


      do 5 it=1,maxit
         if( cvhlm .ge. cvhrm ) then
            if (cvhl .lt. cvhlm) then
               err= cvhl/cvhlm
               ierr=2
               return
            endif
c
            hl=hlm
            cvhl=cvhlm
            hlm=hrm
            hrm=hrm+(hrm-hl)*tau
            cvhlm=cvhrm
            call 
     *  cvrf(hrm,n,x,y,wt,sy,trchrm,diag,din,dout,cvhrm,ierr)  

         else
            if (cvhr .lt. cvhrm) then
               err= cvhrm/cvhr
               ierr=2
               return
            endif
c
            hr=hrm
            cvhr=cvhrm
            hrm=hlm
            hlm=hlm+(hlm-hr)*tau
            cvhrm=cvhlm
            call 
     * cvrf(hlm,n,x,y,wt,sy,trchlm,diag,din,dout,cvhlm,ierr)  
            
         endif
5     continue

**** finished -- take the best h so far
      if( cvhlm .ge. cvhrm) then
         hopt=hrm
         vopt=cvhrm
         tropt=trchrm
      else
         hopt=hlm
         vopt=cvhlm
         tropt=trchlm
      endif
c
10    format(' cv(h) is NOT quasiconvex ',/,
     - ' iter,error,hl,hlm: ',I4,e15.5,2e15.5)
11    format(' cv(h) is NOT quasiconvex ',/,
     - ' iter,error,hr,hrm: ',I4,e15.5,2e15.5)
c
c
      return
      end
