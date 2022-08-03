
c
c
c
c
      subroutine h2d_comb_stab(srcinfo,ndt,targinfo,rdotns,rdotnt,
     1   ndd,dpars,ndz,zpars,ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(8),targinfo(ndt)
      complex *16 zpars(3),u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      zk = zpars(1)
      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = zs*zk*h1*rr

      
      u = zpars(3)*ztmp*rdotns + zpars(2)*zs*h0


      return
      end
c
c
c
c
c
c
      subroutine h2d_sprime_stab(srcinfo,ndt,targinfo,rdotns,rdotnt,
     1   ndd,dpars,ndz,zk,ndi,ipars,u)
      implicit real *8 (a-h,o-z)
      integer ipars
      real *8 dpars,srcinfo(2),targinfo(ndt)
      real *8 rdotns,rdotnt
      complex *16 u,h0,ima,zs,z,zk,h1,gx,gy,ztmp
      data ima/(0.0d0,1.0d0)/
      data zs/(0.0d0,0.25d0)/

      rr2 = (srcinfo(1)-targinfo(1))**2 + (srcinfo(2)-targinfo(2))**2
      rr = dsqrt(rr2)
      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)

      ztmp = -zs*zk*h1*rr

      u = ztmp*rdotnt

      return
      end
c
c
c
c
c
c
      subroutine h2d_transmission_neu_stab(srcinfo,ndt,targinfo,
     1   rdotns,rdotnt,ndd,dpars,ndz,zpars,ndi,ipars,u)
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
      real *8 rdotns,rdotnt
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
      rinv = 1.0d0/rr

      xd = xd*rinv
      yd = yd*rinv


      z = zk*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk

c      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk*zs
c      d2gdxdy = h2*xd*yd*zk*zs
c      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk*zs

      gd0 = h1*rinv*
     1   (srcinfo(7)*targinfo(7) + srcinfo(8)*targinfo(8))
      gd0 = gd0 - h2*rr2*rdotns*rdotnt

c      gd0 = -(d2gdx2*srcinfo(7)*targinfo(7) +
c     1    d2gdxdy*(srcinfo(7)*targinfo(8) + srcinfo(8)*targinfo(7)) + 
c     2    d2gdy2*srcinfo(8)*targinfo(8))

      gs0 = -h1*rr*rdotnt

      u = zk*zs*(zpars(3)*gs0 + zpars(5)*gd0)

      z = zk2*rr
      ifexpon = 1
      call hank103(z,h0,h1,ifexpon)
      h2 = (2/z*h1 - h0)*zk2


      gd1 = h1*rinv*
     1   (srcinfo(7)*targinfo(7) + srcinfo(8)*targinfo(8))
      gd1 = gd1 - h2*rr2*rdotns*rdotnt

c      d2gdx2 = (-h1*rinv + h2*xd*xd)*zk2
c      d2gdxdy = h2*xd*yd*zk2
c      d2gdy2 = (-h1*rinv+h2*yd*yd)*zk2

c      gd1 = -zs*(d2gdx2*srcinfo(7)*targinfo(7) +
c     1    d2gdxdy*(srcinfo(7)*targinfo(8) + srcinfo(8)*targinfo(7)) + 
c     2    d2gdy2*srcinfo(8)*targinfo(8))

c      gx = -zs*zk2*h1*(targinfo(1)-srcinfo(1))*rinv
c      gy = -zs*zk2*h1*(targinfo(2)-srcinfo(2))*rinv
      
      gs1 = gx*targinfo(7) + gy*targinfo(8)
      gs1 = -h1*rr*rdotnt


      u = u+zk2*zs*(zpars(4)*gs1 + zpars(6)*gd1)



      return
      end
c
c
c



