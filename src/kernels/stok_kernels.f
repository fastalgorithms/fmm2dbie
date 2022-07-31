c     We take the following conventions for the Stokes kernels
c
c     For a source y and target x, let r_i = x_i-y_i
c     and let r = sqrt(r_1^2 + r_2^2)
c
c     The Stokeslet, G_{ij}, and its associated pressure tensor, P_j,
c     are
c
c     G_{ij}(x,y) = ((r_i r_j)/(r^2) - delta_{ij}log(r))/(4*pi)
c     P_j(x,y) = (r_j/r^2)/(2*pi)
c
c     The (Type I) stresslet, T_{ijk}, and its associated pressure
c     tensor, PI_{jk}, are
c     
c     T_{ijk}(x,y) = -r_i r_j r_k/ r^4 / pi
c     PI_{jk} = (- delta_{jk}/r^2 + 2 r_j r_k/r^4)/(2*pi)
c
      

      subroutine st2d_slp_vec(nd,srcinfo,ndt,targinfo,ndd,dpars,ndz,
     1     zpars,ndi,ipars,u)
c
c       "single layer" (stokeslet) interaction kernel
c
c       input:
c         srcinfo(2) - double
c           x,y location of source
c         targinfo(2) - double
c           x,y location of target
c         dpars - double
c           dummy parameter
c         zpars - complex
c           dummy parameter
c         ipars(2) - integer
c         i = ipars(1), j = ipars(2). routine returns i,j entry of
c         slp tensor (a 2x2 tensor)
c
c     output:
c         u(4) - double
c     All four entries of G_ij
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars(ndi)
      real *8 dpars,srcinfo(2),targinfo(ndt)
      complex *16 ima,zpars
      data ima/(0.0d0,1.0d0)/

      real *8 dr(2), over4pi, u(nd)
      data over4pi/0.07957747154594767d0/
      data over8pi/0.039788735772973833d0/

      dx=targinfo(1)-srcinfo(1)
      dy=targinfo(2)-srcinfo(2)

      r2=dx**2+dy**2
      rinv2 = (1/r2)*over4pi
      
      diag = -log(r2)*over8pi

      u(1) = diag + dx*dx*rinv2
      u(2) = dx*dy*rinv2
      u(3) = u(2)
      u(4) = diag + dy*dy*rinv2

      return
      end

      subroutine st2d_slp(srcinfo,ndt,targinfo,ndd,dpars,ndz,zpars,ndi,
     1   ipars,val)
c
c       "single layer" (stokeslet) interaction kernel, single entry
c
c       input:
c         srcinfo(2) - double
c           x,y location of source
c         targinfo(2) - double
c           x,y location of target
c         dpars - double
c           dummy parameter
c         zpars - complex
c           dummy parameter
c         ipars(2) - integer
c         i = ipars(1), j = ipars(2). routine returns i,j entry of
c         slp tensor (a 2x2 tensor)
c
c     output:
c         val = G_ij(r)
c     where r is the distance between source and target
c     and i, j are the first two entries of ipars
c     see top of file for def of G_ij
c          
c
      implicit real *8 (a-h,o-z)
      integer ipars(ndi)
      real *8 dpars,srcinfo(2),targinfo(ndt)
      complex *16 ima,zpars
      data ima/(0.0d0,1.0d0)/

      real *8 dr(2), over4pi
      data over4pi/0.07957747154594767d0/      

      i = ipars(1)
      j = ipars(2)
      
      dr(1)=targinfo(1)-srcinfo(1)
      dr(2)=targinfo(2)-srcinfo(2)

      dxi = dr(i)
      dxj = dr(j)

      r2=dr(1)**2+dr(2)**2
      rinv2 = 1/r2
      
      diag = 0
      if (i .eq. j) diag = -log(r2)/2

      val = (diag + dxi*dxj*rinv2)*over4pi

      return
      end




      subroutine st2d_dlp_vec(nd,srcinfo,ndt,targ,ndd,dpars,ndz,zpars,
     1     ndi, ipars,val)
c     f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
      integer ipars(ndi)
      complex *16 :: zpars
      real *8 :: val(nd), over4pi, overpi

      real *8 :: src(2), srcnorm(2)
      data over4pi/0.07957747154594767d0/
      data overpi /0.318309886183790671d0/
c     
c     returns the Stokes double layer kernel
c     
c     D_ij = T_{jik} nu_k
c     
c     where r = |src-targ| and nu is the normal vector at 
c     the source. The output is given ordered by standard
c     linear indexing, ij -> i+(j-1)*2

      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      srcnorm(1)=srcinfo(7)
      srcnorm(2)=srcinfo(8)

c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dprod = dx*srcnorm(1) + dy*srcnorm(2)

      r2=dx**2+dy**2
      rinv4 = -dprod*(1.0d0/r2)**2*overpi

      dxdy = dx*dy*rinv4

      val(1) = dx*dx*rinv4
      val(2) = dxdy
      val(3) = dxdy
      val(4) = dy*dy*rinv4

      return
      end


      subroutine st2d_dlp(srcinfo,ndt,targ,ndd,dpars,ndz, 
     1     zpars,ndi,ipars,val)
c     f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
      integer ipars(ndi)
      complex *16 :: zpars
      real *8 :: val,over4pi,overpi

      real *8 :: src(2), srcnorm(2), dr(2)
      data over4pi/0.07957747154594767d0/
      data overpi /0.318309886183790671d0/
c     
c     returns one entry of the Stokes double layer kernel
c     
c     D_ij = T_jik nu_k 
c     
c     where r = |src-targ| and nu is the normal vector at
c     the source. Returns val=D_ij with i = ipars(1), j = ipars(2)

      i = ipars(1)
      j = ipars(2)
      
      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      srcnorm(1)=srcinfo(7)
      srcnorm(2)=srcinfo(8)

c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dprod = dx*srcnorm(1) + dy*srcnorm(2)

      r2=(dx**2+dy**2)
      rinv4 = -dprod*(1.0d0/r2)**2*overpi

      dr(1) = dx
      dr(2) = dy

      dxi = dr(i)
      dxj = dr(j)
      val = dxi*dxj*rinv4

      return
      end


      subroutine st2d_comb_vec(nd,srcinfo,ndt,targ,ndd,dpars,ndz,zpars,
     1     ndi,ipars,val)
c     f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
      integer ipars(ndi)
      complex *16 :: zpars
      real *8 :: val(nd)

      real *8 :: src(2), srcnorm(2), alpha, beta,over8pi
      data over4pi/0.07957747154594767d0/
      data over8pi/0.039788735772973833d0/      
c     
c     returns the Stokes combined layer kernel
c     
c     K_ij = alpha G_ij + beta D_ij
c     
c     Where alpha = dpars(1), beta = dpars(2), G_ij is the
c     Stokeslet and :
c     
c     D_ij = T_{jik} nu_k
c     
c     where r = |src-targ| and nu is the normal vector at 
c     the source. The output is given ordered by standard
c     linear indexing, ij -> i+(j-1)*2

      alpha = dpars(1)
      beta = dpars(2)
      
      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      srcnorm(1)=srcinfo(7)
      srcnorm(2)=srcinfo(8)

c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dprod = dx*srcnorm(1) + dy*srcnorm(2)

      r2 = dx**2 + dy**2
      rinv2 = 1/r2
      rinv4 = rinv2**2
      diag = -log(r2)*over8pi*alpha
      
      dcomb = (rinv2*alpha-rinv4*beta*dprod*4)*over4pi
      
      dxdy = dx*dy*dcomb

      val(1) = diag + dx*dx*dcomb
      val(2) = dxdy
      val(3) = dxdy
      val(4) = diag + dy*dy*dcomb

      return
      end


      subroutine st2d_comb(srcinfo,ndt,targ,ndd,dpars,ndz, 
     1     zpars,ndi,ipars,val)
c     f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(12), targ(ndt),dpars(ndd)
      integer ipars(ndi)
      complex *16 :: zpars
      real *8 :: val

      real *8 :: src(2), srcnorm(2), dr(2),over4pi
      data over4pi/0.07957747154594767d0/
      data over8pi/0.039788735772973833d0/      
c     
c     returns one entry of the Stokes double layer kernel
c     
c     D_ij = -T_jik n_k = 3*(sum_k (targ_k-src_k) n_k)
c     *(targ_j-src_j)(targ_i-src_i)/r^5
c     
c     where r = |src-targ| and n is the normal vector at
c     the source. Returns val=D_ij with i = ipars(1), j = ipars(2)

      i = ipars(1)
      j = ipars(2)

      alpha = dpars(1)
      beta = dpars(2)
      
      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      srcnorm(1)=srcinfo(7)
      srcnorm(2)=srcinfo(8)

c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dr(1) = dx
      dr(2) = dy
      
      dprod = dx*srcnorm(1) + dy*srcnorm(2)

      r2 = dx**2 + dy**2
      rinv2 = 1/r2
      rinv4 = rinv2**2
      diag = -log(r2)*over8pi*alpha
      
      dcomb = (rinv2*alpha-rinv4*beta*dprod*4)*over4pi
      
      dxi = dr(i)
      dxj = dr(j)
      val = dxi*dxj*dcomb

      if (i .eq. j) val = val + diag

      return
      end



      subroutine st2d_strac_vec(nd,srcinfo,ndt,targinfo,
     1     ndd,dpars,ndz,zpars,ndi,ipars,val)
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(*), targinfo(8),dpars(ndd)
      integer ipars(ndi)
      real *8 :: val(4), targ(2), src(2), targnorm(2),over4pi
      data over4pi/0.07957747154594767d0/
      data overpi /0.318309886183790671d0/
      
c     f2py intent(in) nd,src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val

c     
c     returns the traction of the Stokes single layer kernel
c     
c     t(S)_ij = T_ijk n_k = -(sum_k (targ_k-src_k) n_k)
c     *(targ_j-src_j)(targ_i-src_i)/r^4
c     
c     where r = |src-targ| and n is the normal vector at the
c     target. The output is given ordered by standard
c     linear indexing, ij -> i+(j-1)*2

      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      targ(1)=targinfo(1)
      targ(2)=targinfo(2)
      targnorm(1)=targinfo(7)
      targnorm(2)=targinfo(8)
      
c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dprod = dx*targnorm(1) + dy*targnorm(2)

      r2=dx**2+dy**2
      rinv4 = -dprod*(1.0d0/r2)**2*overpi

      dxdy = dx*dy*rinv4

      val(1) = dx*dx*rinv4
      val(2) = dxdy
      val(3) = dxdy
      val(4) = dy*dy*rinv4

      return
      end


      subroutine st2d_strac(srcinfo,ndt,targinfo,ndd,dpars,
     1     ndz,zpars,ndi,ipars,val)
      implicit real *8 (a-h,o-z)
      real *8 :: srcinfo(*), targinfo(8),dpars(ndd)
      integer ipars(ndi)
      real *8 :: val, targ(2), src(2), targnorm(2), dr(2),over4pi
      complex *16 :: zpars
      data over4pi/0.07957747154594767d0/
      data overpi /0.318309886183790671d0/
      
      
c     f2py intent(in) src,ndt,targ,ndd,dpars,ndz,zpars,ndi,ipars
c     f2py intent(out) val

c     
c     returns one entry of the traction of the Stokes single
c     layer kernel
c     
c     t(S)_ij = T_ijk n_k = -3*(sum_k (targ_k-src_k) n_k)
c     *(targ_j-src_j)(targ_i-src_i)/r^5
c     
c     where r = |src-targ| and n is the normal vector at the
c     target. Returns val=t(S)_ij, with i = ipars(1), j = ipars(2)

      src(1)=srcinfo(1)
      src(2)=srcinfo(2)
      targ(1)=targinfo(1)
      targ(2)=targinfo(2)
      targnorm(1)=targinfo(7)
      targnorm(2)=targinfo(8)
      
c     call prin2('srcinfo=*',srcinfo,12)

      dx=targ(1)-src(1)
      dy=targ(2)-src(2)

      dprod = dx*targnorm(1) + dy*targnorm(2)

      r2=dx**2+dy**2
      rinv4 = -dprod*(1.0d0/r2)**2*overpi

      dr(1) = dx
      dr(2) = dy

      i = ipars(1)
      j = ipars(2)

      dxi = dr(i)
      dxj = dr(j)
      
      val = dxi*dxj*rinv4

      return
      end

