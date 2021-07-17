      
      subroutine l2d_slp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,
     1   ipars,u)
c
c     single layer interaction kernel
c     
c     input:
c      
c     srcinfo(2) - double
c                  x,y location of source
c     targinfo(2) - double
c                  x,y location of target
c     dpars - double
c                  dummy parameter
c     zpars - complex
c                  dummy paramter
c     ipars - integer
c                  dummy parameter
c     
c     output:
c     u = - \log |r| / 2\pi  \, ,
c     where r is the distance between source and target
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(ndt)
      complex *16 zpars
      real *8 u

      u=0

      pi4 = 16*atan(1.0d0)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      
      u = -log(rr2)/pi4

      return
      end


c
c
      subroutine l2d_dlp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,
c
c     double layer interaction kernel
c     
c     input:
c      
c     srcinfo(8) - double
c                   srcinfo(1:2) x,y location of source
c                   srcinfo(7:8) normal vector at source      
c                    
c     targinfo(2) - double
c                  x,y location of target
c     dpars - double
c                  dummy parameter
c     zpars - complex
c                  dummy paramter
c     ipars - integer
c                  dummy parameter
c     
c     output:
c     u = - n(y)\cdot \nabla_y \log |r| / 2\pi  \, ,
c     where r is the distance between source and target
c     \nabla_y is the gradient with respect to the source
c     location, and n(y) is the normal direction at the source
c          
c
      
     1   ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      complex *16 zpars
      real *8 pi,dpars,srcinfo(8),targinfo(2)

      u=0
      thresh=1d-15
      thresh2=thresh*thresh
      
      pi2 = 8*atan(1.0d0)
      
      xdiff = targinfo(1)-srcinfo(1)
      ydiff = targinfo(2)-srcinfo(2)
      rr = xdiff*xdiff+ydiff*ydiff
      if(rr.le. thresh2) goto 1000
      u = (srcinfo(7)*xdiff+srcinfo(8)*ydiff)/rr/pi2
      
 1000 continue
      return
      end
c
c
c
      subroutine l2d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,
     1     ndi,ipars,u)
c
c     normal derivative of single layer interaction kernel
c     
c     input:
c      
c     srcinfo(2) - double
c                   srcinfo(1:2) x,y location of source
c                    
c     targinfo(ndt) - double
c     ndt .ge. 8
c                   targinfo(1:2) x,y location of target
c                   targinfo(7:8) normal direction at target
c     dpars - double
c                  dummy parameter
c     zpars - complex
c                  dummy paramter
c     ipars - integer
c                  dummy parameter
c     
c     output:
c     u = - n(x)\cdot \nabla_x \log |r| / 2\pi  \, ,
c     where r is the distance between source and target
c     \nabla_x is the gradient with respect to the target
c     location, and n(x) is the normal direction at the target
c          
      
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 pi,dpars,srcinfo(8),targinfo(ndt)

      u=0
      thresh=1d-15
      thresh2=thresh*thresh
      
      pi2 = 8*atan(1.0d0)
      
      xdiff = targinfo(1)-srcinfo(1)
      ydiff = targinfo(2)-srcinfo(2)
      rr = xdiff*xdiff+ydiff*ydiff
      if(rr.le. thresh2) goto 1000
      u = -(targinfo(7)*xdiff+targinfo(8)*ydiff)/rr/pi2
      
 1000 continue
      return
      end
c
c
      subroutine l2d_comb(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,
     1     ndi,ipars,u)
c
c     combined layer interaction kernel
c
c     alpha*S + beta*D
c     
c     input:
c      
c     srcinfo(8) - double
c                   srcinfo(1:2) x,y location of source
c                   srcinfo(7:8) normal vector at source      
c                    
c     targinfo(2) - double
c                  x,y location of target
c     dpars(2) - double
c                  alpha=dpars(1),beta=dpars(2)
c     zpars - complex
c                  dummy paramter
c     ipars - integer
c                  dummy parameter
c     
c     output:
c     u = - beta n(y)\cdot \nabla_y \log |r| / 2\pi
c           - alpha \log |r| / 2\pi \, ,
c     where r is the distance between source and target
c     \nabla_y is the gradient with respect to the source
c     location, and n(y) is the normal direction at the source
c     alpha and beta are the coeffs determined by dpars
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars(2),srcinfo(8),targinfo(ndt)
      complex *16 zpars

      alpha=dpars(1)
      beta=dpars(2)
      
      u=0
      thresh=1d-15
      thresh2=thresh*thresh
      
      pi2 = 8*atan(1.0d0)
      pi4 = pi2*2
      
      xdiff = targinfo(1)-srcinfo(1)
      ydiff = targinfo(2)-srcinfo(2)
      rr = xdiff*xdiff+ydiff*ydiff
      if(rr.le. thresh2) goto 1000
      u = beta*(srcinfo(7)*xdiff+srcinfo(8)*ydiff)/rr/pi2
      u = u-alpha*log(rr)/pi4
      

 1000 continue

      return
      end
c
