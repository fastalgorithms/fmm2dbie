
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

      ztmp = -zs*zk*h1*rr

      
      u = zpars(3)*ztmp*rdotns + zpars(2)*zs*h0


      return
      end
c
c
