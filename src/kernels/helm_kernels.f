      subroutine h2d_slp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,ndi,
     1   ipars,u)
c
c       single layer interaction kernel
c
c       input:
c         srcinfo(2) - double
c           x,y location of source
c         targinfo(2) - double
c           x,y location of target
c         dpars - double
c           dummy parameter
c         zk - complex
c           Helmholtz parameter
c         ipars - integer
c           dummy parameter
c
c       output:
c         u = i/4 H_{0}(zk*|r|) \, ,
c         where r is the distance between source and target
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(ndt)
      complex *16 u,h0,ima,zs,z,zk,h1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1

      call hank103(z,h0,h1,ifexpon)
      u = zs*h0
      

      return
      end


c
c
c
c
c
c
      subroutine h2d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,
     1   ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(ndt)
      complex *16 u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(targinfo(1)-srcinfo(1))
      gy = ztmp*(targinfo(2)-srcinfo(2))
      
      u = gx*targinfo(7) + gy*targinfo(8)

      return
      end
c
c
c
c
c
      subroutine h2d_dlp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,
     1   ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(8),targinfo(2)
      complex *16 u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = gx*srcinfo(7) + gy*srcinfo(8)

      return
      end
c
c
c
c
c
      subroutine h2d_comb(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,
     1   ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(8),targinfo(ndt)
      complex *16 zpars(3),u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = zpars(3)*(gx*srcinfo(7) + gy*srcinfo(8)) + zpars(2)*zs*h0


      return
      end
c
c
c        transmission kernels
c
c
      subroutine h2d_transmission_dir(srcinfo,ndt,targinfo,ndd,dpars,
     1   ndz,zpars,ndi,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1} + beta S_{k2} + gamma D_{k1} + delta D_{k2}
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(8),targinfo(2)
      complex *16 zpars(6),u,h0,ima,zs,z,zk,zk2,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)

      rinv = 1/rr


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = zpars(5)*(gx*srcinfo(7) + gy*srcinfo(8)) + zpars(3)*zs*h0

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk2*h1*rinv

      gx = ztmp*(srcinfo(1)-targinfo(1))
      gy = ztmp*(srcinfo(2)-targinfo(2))
      
      u = u+zpars(6)*(gx*srcinfo(7) + gy*srcinfo(8)) + zpars(4)*zs*h0


      return
      end
c
c
c
c
c
      subroutine h2d_transmission_neu(srcinfo,ndt,targinfo,ndd,dpars,
     1   ndz,zpars,ndi,ipars,u)
c
c
c         The kernel of interaction is given by
c           alpha S_{k1}' + beta S_{k2}' + gamma D_{k1}' + delta D_{k2}'
c         
c          zpars(1) = k1
c          zpars(2) = k2
c          zpars(3:6) = alpha,beta,gamma,delta
c          
c

      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(8),targinfo(8)
      complex *16 zpars(6),u,h0,ima,zs,z,zk,h1,gx,gy,h2,zk2
      complex *16 d2gdx2,d2gdy2,d2gdxdy,ztmp
      complex *16 gd0,gs0,gd1,gs1
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      zk2 = zpars(2)

      xd = targinfo(1) - srcinfo(1)
      yd = targinfo(2) - srcinfo(2)
      
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      rinv = 1/rr

      xd = xd*rinv
      yd = yd*rinv


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk

      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk*zs
      d2gdxdy = h2*xd*yd*zk*zs
      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk*zs

      gd0 = -(d2gdx2*srcinfo(7)*targinfo(7) +
     1    d2gdxdy*(srcinfo(7)*targinfo(8) + srcinfo(8)*targinfo(7)) + 
     2    d2gdy2*srcinfo(8)*targinfo(8))

      gx = -zs*zk*h1*(targinfo(1)-srcinfo(1))*rinv
      gy = -zs*zk*h1*(targinfo(2)-srcinfo(2))*rinv

      gs0 = gx*targinfo(7) + gy*targinfo(8)


      u = zpars(3)*gs0 + zpars(5)*gd0

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk2

      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk2
      d2gdxdy = h2*xd*yd*zk2
      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk2

      gd1 = -zs*(d2gdx2*srcinfo(7)*targinfo(7) +
     1    d2gdxdy*(srcinfo(7)*targinfo(8) + srcinfo(8)*targinfo(7)) + 
     2    d2gdy2*srcinfo(8)*targinfo(8))

      gx = -zs*zk2*h1*(targinfo(1)-srcinfo(1))*rinv
      gy = -zs*zk2*h1*(targinfo(2)-srcinfo(2))*rinv
      
      gs1 = gx*targinfo(7) + gy*targinfo(8)


      u = u+zpars(4)*gs1 + zpars(6)*gd1



      return
      end
c
c
c
      subroutine l2d_dlp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zk,
     1   ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 pi,dpars,srcinfo(8),targinfo(2)
      complex *16 u,ima,z,zinv
      complex *16 dipstr
      data ima/(0.0d0,1.0d0)/

      dipstr = dcmplx(srcinfo(7),srcinfo(8))
      pi = 4.0d0*datan(1.0d0)
c
      xdiff = targinfo(1)-srcinfo(1)
      ydiff = targinfo(2)-srcinfo(2)
      z = dcmplx(xdiff,ydiff)
      if(abs(z).le.1.0d-16) goto 1000
      zinv = 1/z
      u= zinv*dipstr/(2*pi)
c
1000  continue
      return
      end
c
c

