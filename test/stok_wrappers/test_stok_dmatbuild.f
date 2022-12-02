c
c     unit test of the complex matrix builder (scalar problem)
c      
      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      integer, allocatable :: adjs(:,:), ipiv(:)
      real *8 ts(2)
      character *100 fname
      integer ipars(2)
      real *8, allocatable :: ubdry(:,:), uin(:,:), qwts(:)
      real *8, allocatable :: potbdry(:,:), potin(:,:), potin2(:,:)
      real *8, allocatable :: sysmat(:,:)
      
      complex *16 zk, ima

      integer, allocatable :: norders(:),ixys(:),iptype(:)
      integer, allocatable :: ixyso(:),nfars(:)

      integer, allocatable :: ich_id(:),inode_id(:)
      real *8, allocatable :: ts_targ(:)
      real *8, allocatable :: ab(:,:)
      real *8 :: errs(500), rres(500)
      real *8 xyz_out(2),xyz_in(2),gmat(2,2),mu(2)
      real *8, allocatable :: sigma(:,:), sigma2(:,:)
      complex * 16 zpars
      real *8 dpars(2), alpha, beta
      integer opdims(2), info
      external fstarn_simple
      procedure (), pointer :: fker
      external st2d_slp_wrap, st2d_dlp_wrap, st2d_comb_wrap
      data ima /(0d0,1d0)/

      real *8 ucc(2), ucc2(2)
      
      call prini(6,13)

      done = 1
      pi = atan(done)*4

      alpha=1
      beta=1

      
      nchmax = 20000
      k = 16
      npts_max = nchmax*k

      allocate(srcvals(8,npts_max),srccoefs(6,npts_max),ab(2,nchmax))
      allocate(norders(nchmax),ixys(nchmax+1),iptype(nchmax))
      allocate(adjs(2,nchmax))

      irefinel = 0
      irefiner = 0
      rlmax = huge(1d0)
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
      eps = 1.0d-6
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1     ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi,ipars,nover,
     2     k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,
     3     adjs,ier)

      call prinf("nch *",nch,1)
      call prinf("npts *",npts,1)      


      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0
      mu(1) = 1.0d0
      mu(2) = 1.1d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0

      ntin = 1
      
      allocate(sigma(2,npts),ubdry(2,npts),uin(2,ntin))
      allocate(sigma2(2,npts))

c     get boundary data
      
      do i=1,npts
        call st2d_slp_wrap(xyz_out,2,srcvals(1,i),0,dpars,1,zpars,0,
     1        ipars,gmat)
        ubdry(1,i) = gmat(1,1)*mu(1) + gmat(1,2)*mu(2)
        ubdry(2,i) = gmat(2,1)*mu(1) + gmat(2,2)*mu(2)
      enddo

c     test data
      
      call st2d_slp_wrap(xyz_out,2,xyz_in,0,dpars,1,zpars,0,
     1     ipars,gmat)

      uin(1,1) = gmat(1,1)*mu(1)+gmat(1,2)*mu(2)
      uin(2,1) = gmat(2,1)*mu(1)+gmat(2,2)*mu(2)
      

c     build system matrix (1/2 beta*I + alpha*S + beta*D + W)

      allocate(sysmat(2*npts,2*npts),ipiv(2*npts))

      ising = 0
      iquadtype = 0
      ifrobust = 0
      opdims(1) = 2
      opdims(2) = 2
      ndz = 0
      ndd = 2
      ndi = 0
      dpars(1) = alpha
      dpars(2) = beta
      fker => st2d_comb_wrap
      call cpu_time(t1)
      call dchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquadtype,ifrobust,
     3     sysmat,ier)
      call cpu_time(t2)

      
      write(*,*) 'time matbuild', t2-t1
      write(*,*) 'npts ', npts

      allocate(qwts(npts))
      call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,qwts)      
      
      do i = 1,npts
         ii = 2*(i-1)+1
         sysmat(ii,ii) = sysmat(ii,ii)+0.5d0*beta
         ii = ii + 1
         sysmat(ii,ii) = sysmat(ii,ii)+0.5d0*beta
         sigma2(1,i) = ubdry(1,i)
         sigma2(2,i) = ubdry(2,i)
         do j = 1,npts
            ii = 2*(i-1)+1
            jj = 2*(j-1)+1
            sysmat(jj,ii) = sysmat(jj,ii) + qwts(i)
            ii = ii+1
            jj = jj+1
            sysmat(jj,ii) = sysmat(jj,ii) + qwts(i)
         enddo
      enddo

      call cpu_time(t1)
      nrhs=1
      npts2=2*npts
      call dgesv(npts2, nrhs, sysmat, npts2, ipiv, sigma2, npts2,
     1     info)
      call cpu_time(t2)

      ucc2(:)=0
      do i = 1,npts
         
         ucc2(1) = ucc2(1) + qwts(i)*sigma2(1,i)
         ucc2(2) = ucc2(2) + qwts(i)*sigma2(2,i)
      enddo

      
      call prinf('info *',info,1)
      call prin2('sigma2 *',sigma2,10)

      write(*,*) 'time solve with LAPACK', t2-t1
      write(*,*) 'info ', info

      dpars(1) = alpha
      dpars(2) = beta
      numit=100
      ifinout=0
      ifnocorr=0
      eps_gmres=1d-14
      call cpu_time(t1)
      call stok_comb_dir_solver_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,eps,dpars,
     2     ifnocorr,numit,
     2     ifinout,ubdry,eps_gmres,niter,errs,rres,sigma,ucc)
      call cpu_time(t2)

      call prin2('sigma *',sigma,10)

      
      write(*,*) 'time solve with GMRES', t2-t1
      

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

      dpars(1) = alpha
      dpars(2) = beta

      allocate(potbdry(2,npts),potin(2,ntin),potin2(2,ntin))
      
      ndtarg2 = 2
      n1 = 1
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg2,n1,xyz_in,ich_id,
     2     ts_targ,eps,dpars,sigma,potin)
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg2,n1,xyz_in,ich_id,
     2     ts_targ,eps,dpars,sigma2,potin2)


      potin(1,1) = potin(1,1)+ucc(1)
      potin(2,1) = potin(2,1)+ucc(2)
      potin2(1,1) = potin2(1,1)+ucc2(1)
      potin2(2,1) = potin2(2,1)+ucc2(2)

      call prin2('pot in *',potin,2)
      call prin2('pot in 2 *',potin2,2)
      call prin2('u in *',uin,2)

      write(*,*) abs(potin(1,1)-uin(1,1))
      write(*,*) abs(potin(2,1)-uin(2,1))
      write(*,*) abs(potin2(1,1)-uin(1,1))
      write(*,*) abs(potin2(2,1)-uin(2,1))

      nsuccess = 0
      if (abs(potin(1,1)-uin(1,1)) .lt. eps) nsuccess=1
      ntests=1
      
c      open(unit=33,file='../../print_testres.txt',access='append')
c      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
c     1  ' out of ',ntests,' in helm_wrappers matrix build testing suite'
c      close(33)

      stop
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



      subroutine st2d_dlp_wrap(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(8), y(ndtarg), dpars(ndd)
      complex *16 zpars
      integer ipars
      real *8 f(2,2)

      call st2d_dlp_vec(4,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      return
      end

      subroutine st2d_slp_wrap(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(8), y(ndtarg), dpars(ndd)
      complex *16 zpars
      integer ipars
      real *8 f(2,2)

      call st2d_slp_vec(4,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      return
      end


      subroutine st2d_comb_wrap(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(8), y(ndtarg), dpars(ndd)
      complex *16 zpars
      integer ipars
      real *8 f(2,2)

      call st2d_comb_vec(4,x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,
     1     ipars,f)
      return
      end

      
