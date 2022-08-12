      implicit real *8 (a-h,o-z)
      real *8 srcinfo(8),targinfo(8)
      complex *16 z1,z2,z3,z4,zpars(6),ima,zk
      data ima/(0.0d0,1.0d0)/
      call prini(6,13)

      do i=1,8
        srcinfo(i) = hkrand(0)
        targinfo(i) = hkrand(0)
      enddo
      do i=1,6
        zpars(i) = hkrand(0) + ima*hkrand(0)
      enddo

      x = targinfo(1)-srcinfo(1)
      y = targinfo(2)-srcinfo(2)
      r2 = x**2 + y**2
      rdotns = (x*srcinfo(7) + y*srcinfo(8))/r2
      rdotnt = (x*targinfo(7) + y*targinfo(8))/r2
c 
c  test helmholtz combined field rep
c
      ndt = 8
      ndz = 3
      call h2d_comb(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,
     1 ipars,z1)
      
      call h2d_comb_stab(srcinfo,ndt,targinfo,rdotns,rdotnt,
     1 ndd,dpars,ndz,zpars,ndi,ipars,z2)
      call prin2_long('z1=*',z1,2)
      call prin2_long('z2=*',z2,2)
      err1 = abs(z1-z2)/abs(z1)
      call prin2('error in h2d_comb_stab=*',err1,1)

      ndz = 1
      call h2d_sprime(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,
     1 ipars,z1)
      
      call h2d_sprime_stab(srcinfo,ndt,targinfo,rdotns,rdotnt,
     1 ndd,dpars,ndz,zpars,ndi,ipars,z2)
      call prin2_long('z1=*',z1,2)
      call prin2_long('z2=*',z2,2)
      err1 = abs(z1-z2)/abs(z1)
      call prin2('error in h2d_sprime_stab=*',err1,1)

      ndz = 6 
      call h2d_transmission_neu(srcinfo,ndt,targinfo,ndd,dpars,
     1  ndz,zpars,ndi,ipars,z1)
      
      call h2d_transmission_neu_stab(srcinfo,ndt,targinfo,rdotns,rdotnt,
     1 ndd,dpars,ndz,zpars,ndi,ipars,z2)
      call prin2_long('z1=*',z1,2)
      call prin2_long('z2=*',z2,2)
      err1 = abs(z1-z2)/abs(z1)
      call prin2('error in h2d_transmission_neu_stab=*',err1,1)


      stop
      end
      
