
      subroutine test_lege_extras()
      implicit real *8 (a-h,o-z)


      real *8 pols(100), x(100), w(100), u(10000), v(10000)
      complex *16 zpols(100), z, ima
      data ima / (0d0,1d0) /

      call prini(6,13)

      x1 = 0.75d0
      z = x1

      n = 16
      call legepols(x1,n,pols)
      call legepolz(z,n,zpols)

      call prin2('pols *',pols,n+1)
      call prin2('zpols *',zpols,2*n+2)    

      th = 3d0/8
      rho = 1.5d0
      z = rho*exp(ima*th)
      z = (z+1/z)/2
      call legepolz(z,n,zpols)      

      call prin2('rho *',rho,1)
      call prin2('z *',z,2)
      call prin2('zpols *',zpols,2*n+2)    


      itype = 1
      k=64
      call legeexps(itype,k,x,u,v,w)
      wmax = 0
      do i = 1,k
         wmax = max(w(i),wmax)
      enddo

      call prin2('wmax *',wmax,1)
      eps = 1d-12

      rho = (wmax/(eps*2.5))**(1d0/k)

      z = rho*exp(ima*th)
      z = (z+1/z)/2

      call legepolz(z,n,zpols)      

      call prin2('rho *',rho,1)
      call prin2('z *',z,2)
      call prin2('zpols *',zpols,2*n+2)    
      
      stop
      end

