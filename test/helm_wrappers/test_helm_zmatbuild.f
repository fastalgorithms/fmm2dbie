c
c     unit test of the complex matrix builder (scalar problem)
c      
      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      integer, allocatable :: adjs(:,:), ipiv(:)
      real *8 errs(6),ts(2)
      character *100 fname
      integer ipars(2)
      complex *16, allocatable :: ubdry(:), uin(:)
      complex *16, allocatable :: potbdry(:), potin(:)
      complex *16, allocatable :: sysmat(:,:)
      
      complex *16 zk

      integer, allocatable :: norders(:),ixys(:),iptype(:)
      integer, allocatable :: ixyso(:),nfars(:)

      integer, allocatable :: ich_id(:),inode_id(:)
      real *8, allocatable :: ts_targ(:)
      real *8, allocatable :: ab(:,:)
      real *8 xyz_out(2),xyz_in(2)
      complex *16, allocatable :: sigma(:)
      complex * 16 zpars(3)
      real *8 dpars(2)
      integer opdims(2), info
      external fstarn_simple
      procedure (), pointer :: fker
      external h2d_slp, h2d_dlp, h2d_comb


      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nchmax = 20000
      k = 16
      npts_max = nchmax*k

      allocate(srcvals(8,npts_max),srccoefs(6,npts_max),ab(2,nchmax))
      allocate(norders(nchmax),ixys(nchmax+1),iptype(nchmax))
      allocate(adjs(2,nchmax))

      irefinel = 0
      irefiner = 0
      rlmax = 1.0d10
      ta = -pi
      tb = pi
      ifclosed = 1
      rlmaxe = 0

      ndd = 2
      ndi = 1
      ndz = 0

      dpars(1) = 1
      dpars(2) = 0.25d0
      ipars(1) = 5
      
      nover = 1
      nch = 0
      ier = 0
      eps = 1.0d-8
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1  ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi,ipars,nover,
     2  k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,adjs,
     3     ier)

      call prinf("nch *",nch,1)
      call prinf("npts *",npts,1)      



      zk = 1.0d0
      zpars(1) = zk
      zpars(2) = 1.0d0
      zpars(3) = 0.0d0

      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0

      ntin = 1
      
      allocate(sigma(npts),ubdry(npts),uin(ntin))

c     get boundary data
      
      do i=1,npts
        call h2d_slp(xyz_out,2,srcvals(1,i),0,dpars,1,zpars,0,
     1     ipars,ubdry(i))
      enddo

c     test data
      
      call h2d_slp(xyz_out,2,xyz_in,0,dpars,1,zpars,0,
     1     ipars,uin)


c     build system matrix (-1/2 I + D)

      allocate(sysmat(npts,npts),ipiv(npts))

      ising = 0
      iquadtype = 0
      ifrobust = 0
      opdims(1) = 1
      opdims(2) = 1
      ndz = 1
      ndd = 0
      ndi = 0
      zpars(1) = zk
      fker => h2d_dlp
      call cpu_time(t1)
      call zchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquadtype,ifrobust,
     3     sysmat,ier)
      call cpu_time(t2)

      
      write(*,*) 'time matbuild', t2-t1
      write(*,*) 'npts ', npts

      
      do i = 1,npts
         sysmat(i,i) = sysmat(i,i)-0.5d0
         sigma(i) = ubdry(i)
      enddo

      call cpu_time(t1)
      nrhs=1
      call zgesv(npts, nrhs, sysmat, npts, ipiv, sigma, npts,
     1     info)
      call cpu_time(t2)

      write(*,*) 'time solve with LAPACK', t2-t1
      write(*,*) 'info ', info
      

      ndtarg = 2
      allocate(targs(ndtarg,npts))
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
      enddo

      allocate(ich_id(npts),ts_targ(npts))
      do i=1,npts
        ich_id(i) = -1
        ts_targ(i) = 0
      enddo

      zpars(1) = zk
      zpars(2) = 0.0d0
      zpars(3) = 1.0d0

      allocate(potbdry(npts),potin(ntin))
      
      call lpcomp_helm_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ich_id,
     2     ts_targ,eps,zpars,sigma,potbdry)

      ndtarg2 = 2
      n1 = 1
      call lpcomp_helm_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg2,n1,xyz_in,ich_id,
     2     ts_targ,eps,zpars,sigma,potin)

      call prin2('pot in *',potin,2)
      call prin2('u in *',uin,2)

      write(*,*) abs(potin(1)-uin(1))

      nsuccess = 0
      if (abs(potin(1)-uin(1)) .lt. eps) nsuccess=1
      ntests=1
      
      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in helm_wrappers matrix build testing suite'
      close(33)

      return      
      end



c
c
c
c
c
      subroutine fstarn_simple(t,ndd,dpars,ndz,zpars,ndi,ipars,x,
     1  y,dxdt,dydt,dxdt2,dydt2)
c
c   curve parametrization for simple star shaped domain
c
c       r(t) = r1 + r2*cos(n*t)
c       x(t) = r(t) cos(t)
c       y(t) = r(t) sin(t)
c   
c   where r1 = dpars(1), r2 = dpars(2) and n = ipars(1)
c     
c


      implicit real *8 (a-h,o-z)
      real *8 dpars(2)
      integer ipars
      complex *16 zpars,z,zd,zd2,zexp,ima,ze1,ze2,ze3
      data ima/(0.0d0,1.0d0)/

      rsc = 1.0d0

      n = ipars
      r1 = dpars(1)*rsc
      r2 = dpars(2)*rsc

      r = r1 + r2*cos(n*t)


      drdt = -r2*(n+0.0d0)*sin(n*t)
      d2rdt2 = -r2*(n+0.0d0)*(n+0.0d0)*cos(n*t)

      x = r*cos(t)
      y = r*sin(t)
      
      dxdt = drdt*cos(t) - r*sin(t)
      dydt = drdt*sin(t) + r*cos(t)
      
      dxdt2 = d2rdt2*cos(t) - 2*drdt*sin(t) - r*cos(t) 
      dydt2 = d2rdt2*sin(t) + 2*drdt*cos(t) - r*sin(t) 
      
      z = (r1 + r2*cos((n+0.0d0)*t))*exp(ima*t)
      zd = -(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t) + ima*z
      zd2 = ima*zd - (n+0.0d0)**2*r2*cos((n+0.0d0)*t)*exp(ima*t) -
     1    ima*(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t)
      ze1 = exp(ima*t)
      ze2 = exp(ima*t*(n+1))
      ze3 = exp(ima*t*(-n+1))
      z = r1*ze1 + 0.5d0*r2*(ze2+ze3)
      zd = ima*(r1*ze1 + 0.5d0*r2*(n+1)*ze2 + 0.5d0*r2*(-n+1)*ze3)
      zd2 = -(r1*ze1 + 0.5d0*r2*(n+1)**2*ze2 + 0.5d0*r2*(-n+1)**2*ze3)
      
      x = real(z)
      y = imag(z)
      dxdt = real(zd)
      dydt = imag(zd)
      dxdt2 = real(zd2)
      dydt2 = imag(zd2)


      return
      end



