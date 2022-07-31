      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      integer, allocatable :: adjs(:,:)
      real *8 errs(6),ts(2)
      character *100 fname
      real *8, allocatable :: uval(:,:),tracval(:,:)
      real *8, allocatable :: pot(:,:),potslp(:,:),potdlp(:,:)
      real *8, allocatable :: pot2(:,:),potslp2(:,:),potdlp2(:,:)
      real *8, allocatable :: slpmat(:,:), dlpmat(:,:)

      integer opdims(2)
      
      integer, allocatable :: norders(:),ixys(:),iptype(:)

      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_targ(:)
      real *8, allocatable :: ab(:,:)
      real *8 xyz_out(2),xyz_in(2)
      integer ipars(2)
      complex *16 zpars
      real *8 dpars(2), gmat(2,2), tracmat(2,2), musrc(2)

      real *8 potin(2), potintrue(2), potslpin(2), potdlpin(2)
      real *8 ts_targ_in(1)
      integer ich_id_in(1)
      

      external fstarn_simple, st2d_slp_wrap, st2d_dlp_wrap



      call prini(6,13)


      call stok_kernel_test(ipass)

      write(*,*) ipass

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
      
      dpars(1) = 1d0
      dpars(2) = 0.2d0

      ipars(1) = 3
      nover = 3
      nch = 0
      ier = 0
      eps = 1.0d-12
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1  ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi,ipars,nover,
     2  k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,adjs,
     3  ier)

      eps = 1d-8

      xyz_out(1) = 3.17d0
      xyz_out(2) = 1.3d0

      xyz_in(1) = 1.199999d0
      write(*,*) 1.2d0-xyz_in(1)
      xyz_in(2) = 0

      allocate(targs(2,npts))
      allocate(wts(npts))
      call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,wts)


      allocate(pot(2,npts),potslp(2,npts),potdlp(2,npts))
      allocate(uval(2,npts),tracval(2,npts))
      allocate(ich_id(npts),ts_targ(npts))

      musrc(1) = 1.1d0
      musrc(2) = -.23d0
      
      do i=1,npts
         call st2d_slp_vec(4,xyz_out,2,srcvals(1,i),
     1        0,dpars,0,zpars,0,ipars,gmat)
         uval(1,i) = gmat(1,1)*musrc(1)+gmat(1,2)*musrc(2)
         uval(2,i) = gmat(2,1)*musrc(1)+gmat(2,2)*musrc(2)
         call st2d_strac_vec(4,xyz_out,8,srcvals(1,i),
     1        0,dpars,0,zpars,0,ipars,tracmat)
         tracval(1,i) = tracmat(1,1)*musrc(1)+tracmat(1,2)*musrc(2)
         tracval(2,i) = tracmat(2,1)*musrc(1)+tracmat(2,2)*musrc(2)
      enddo

      dmaxu1 = 0
      dmaxu2 = 0
      do i = 1,npts
         dmaxu1 = max(dmaxu1,abs(uval(1,i)))
         dmaxu2 = max(dmaxu2,abs(uval(2,i)))
      enddo
      write(*,*) 'max uval', dmaxu1, dmaxu2

c     test just at interior target
      
      ndtarg = 2
      do i=1,1
        ich_id_in(i) = -1
        ts_targ_in(i) = 0
      enddo

      dpars(1) = 1
      dpars(2) = 0

      ntarg1 = 1

      call cpu_time(t1)
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg1,xyz_in,
     2     ich_id_in,ts_targ_in,eps,dpars,tracval,potslpin)

      dpars(1) = 0
      dpars(2) = 1
      
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg1,xyz_in,
     2     ich_id_in,ts_targ_in,eps,dpars,uval,potdlpin)

      call cpu_time(t2)
      tlpcomp = t2-t1

      potin(1) = potslpin(1)+potdlpin(1)
      potin(2) = potslpin(2)+potdlpin(2)

      call st2d_slp_vec(4,xyz_out,2,xyz_in,
     1     0,dpars,0,zpars,0,ipars,gmat)
      potintrue(1) = gmat(1,1)*musrc(1)+gmat(1,2)*musrc(2)
      potintrue(2) = gmat(2,1)*musrc(1)+gmat(2,2)*musrc(2)


      errl2 = (potintrue(1)-potin(1))**2 + (potintrue(2)-potin(2))**2
      rell2 =  potintrue(1)**2 + potintrue(2)**2

      err = sqrt(errl2/rell2)
      
      call prin2('error in greens identity, interior pt=*',err,1)
      
      
c     test on boundary
      

      i1 = 0
      if(err.lt.1.0d-2) i1 = 1
      ndtarg = 2
      
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
      enddo

      do i=1,npts
        ich_id(i) = -1
        ts_targ(i) = 0
      enddo

      call get_chunk_id_ts(nch,norders,ixys,iptype,npts, 
     1         ich_id,ts_targ)


      dpars(1) = 1
      dpars(2) = 0
      
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2     ich_id,ts_targ,eps,dpars,tracval,potslp)

      dpars(1) = 0
      dpars(2) = 1
      
      call lpcomp_stok_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2     ich_id,ts_targ,eps,dpars,uval,potdlp)

      call cpu_time(t2)
      tlpcomp = t2-t1

c
c      compute error
c
      errl2 = 0
      rl2 = 0
      do i=1,npts
        pot(1,i) = (potslp(1,i) + potdlp(1,i))*2
        pot(2,i) = (potslp(2,i) + potdlp(2,i))*2
        errl2 = errl2 + abs(uval(1,i)-pot(1,i))**2*wts(i)
        errl2 = errl2 + abs(uval(2,i)-pot(2,i))**2*wts(i)
        rl2 = rl2 + abs(uval(1,i))**2*wts(i)
        rl2 = rl2 + abs(uval(2,i))**2*wts(i)

        if(i.lt.5) then
           write(*,*) potslp(1,i), potdlp(1,i)
c           write(*,*) i
c           write(*,*) uval(1,i)/pot(1,i), uval(1,i), pot(1,i)
c           write(*,*) uval(2,i)/pot(2,i), uval(2,i), pot(2,i)
        endif
      enddo


      err = sqrt(errl2/rl2)

      call prin2('error in greens identity, on bdry=*',err,1)

      i1 = 0
      if(err.lt.1.0d-2) i1 = 1

      opdims(1)=2
      opdims(2)=2
      ising = 0
      iquad = 1
      ifrobust = 0
      npts2 = 2*npts
      allocate(slpmat(npts2,npts2),dlpmat(npts2,npts2))
      
      call dchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,st2d_slp_wrap,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquad,ifrobust,
     3     slpmat,ier)

      call dchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,st2d_dlp_wrap,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquad,ifrobust,
     3     dlpmat,ier)

      allocate(potdlp2(2,npts),potslp2(2,npts),pot2(2,npts))
      
      call test_matvec(npts2,npts2,dlpmat,uval,potdlp2)
      call test_matvec(npts2,npts2,slpmat,tracval,potslp2)


      errl2 = 0
      rl2 = 0
      do i=1,npts
         pot2(1,i) = (potslp2(1,i) + potdlp2(1,i))*2
         pot2(2,i) = (potslp2(2,i) + potdlp2(2,i))*2
      enddo

      call test_wl2err(2,npts,wts,uval,pot2,errl2,rell2,relerr)

      call prin2('rel error in greens id, on bdry (ggq)=*',relerr,1)

      call test_wl2err(2,npts,wts,potslp2,potslp,errl2,rell2,relerr)

      call prin2('l2 diff S[trac] (ggq vs fmm)=*',relerr,1)
      
      call test_wl2err(2,npts,wts,potdlp2,potdlp,errl2,rell2,relerr)

      call prin2('l2 diff D[u] (ggq vs fmm)=*',relerr,1)
      
      
      nsuccess = i1
      ntests = 1

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in helm_wrappers testing suite'
      close(33)

      stop
      end


      subroutine test_l2err(nd,n,xt,x,errl2,rell2,relerr)
      implicit real *8 (a-h,o-z)
      real *8 xt(nd,n), x(nd,n)

      errl2 = 0
      rell2 = 0

      do i = 1,n
         do j = 1,nd
            errl2 = errl2 + (xt(j,i)-x(j,i))**2
            rell2 = rell2 + (xt(j,i))**2
         enddo
      enddo

      errl2 = sqrt(errl2)
      rell2 = sqrt(rell2)
      relerr = errl2/rell2
      
      return
      end

      subroutine test_wl2err(nd,n,w,xt,x,errl2,rell2,relerr)
      implicit real *8 (a-h,o-z)
      real *8 xt(nd,n), x(nd,n), w(n)

      errl2 = 0
      rell2 = 0

      do i = 1,n
         do j = 1,nd
            errl2 = errl2 + (xt(j,i)-x(j,i))**2*w(i)
            rell2 = rell2 + (xt(j,i))**2*w(i)
         enddo
      enddo

      errl2 = sqrt(errl2)
      rell2 = sqrt(rell2)
      relerr = errl2/rell2
      
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

      subroutine test_matvec(m,n,a,x,y)
      implicit real *8 (a-h,o-z)
      real *8 a(m,n), x(n), y(m)

      do i = 1,m
         y(i)=0
      enddo
      
      do j = 1,n
         xj = x(j)
         do i = 1,m
            y(i) = y(i) + a(i,j)*xj
         enddo
      enddo

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



      subroutine stok_kernel_test(ipass)
      implicit real *8 (a-h,o-z)
      real *8 srcinfo(8), targinfo(8), gmat(2,2), tracmat(2,2)
      real *8 u(2), mu(2), nu(2), trac(2), trac2(2)
      real *8 pot(2), pre, grad(2,2), e(2,2)
      real *8 dpars(2), stoklet(2), strslet(2)
      integer ipars(2)

      ipass=1

      pi = 4*atan(1.0d0)
      
      srcinfo(1) = 0.2d0
      srcinfo(2) = 0.3d0


      srcinfo(7) = cos(0.1d0)
      srcinfo(8) = sin(0.1d0)

      targinfo(1) = -0.15d0
      targinfo(2) = 0.2d0

      targinfo(7) = cos(0.2d0)
      targinfo(8) = sin(0.2d0)

      mu(1) = 0.3d0
      mu(2) = -0.1d0

c
c     SLP tests
c
      
      call st2d_slp_vec(4,srcinfo,2,targinfo,0,dpars,0,zpars,
     1     0,ipars,gmat)

      
      u(1) = gmat(1,1)*mu(1) + gmat(1,2)*mu(2)
      u(2) = gmat(2,1)*mu(1) + gmat(2,2)*mu(2)

      thresh = 1d-15
      nd1 = 1
      ns1 = 1
      nt1 = 1
      pot(1)=0
      pot(2)=0
      pre = 0
      grad(1,1) = 0
      grad(2,1) = 0
      grad(1,2) = 0
      grad(2,2) = 0
      call st2ddirectstokg(nd1,srcinfo,mu,ns1,targinfo,nt1,
     1     pot,pre,grad,thresh)
      pot(1) = pot(1)/(2*pi)
      pot(2) = pot(2)/(2*pi)

      do j = 1,2
         do i = 1,2
            ipars(1) = i
            ipars(2) = j
            call st2d_slp(srcinfo,2,targinfo,0,dpars,0,zpars,
     1           2,ipars,g1)
            if (abs(g1-gmat(i,j)) .gt. 1d-14) ipass=0
         enddo
      enddo

      if (sqrt( (pot(1)-u(1))**2 + (pot(2)-u(2))**2) .gt. 1d-14)
     1     ipass=0
      
      
c
c     DLP tests
c
      
      call st2d_dlp_vec(4,srcinfo,2,targinfo,0,dpars,0,zpars,
     1     0,ipars,gmat)

      u(1) = gmat(1,1)*mu(1) + gmat(1,2)*mu(2)
      u(2) = gmat(2,1)*mu(1) + gmat(2,2)*mu(2)

      thresh = 1d-15
      nd1 = 1
      ns1 = 1
      nt1 = 1
      nu(1) = srcinfo(7)
      nu(2) = srcinfo(8)
      ifstoklet = 0
      istress=1
      pot(1)=0
      pot(2)=0
      pre = 0
      grad(1,1) = 0
      grad(2,1) = 0
      grad(1,2) = 0
      grad(2,2) = 0
      call st2ddirectstokstrsg(nd1,srcinfo,ifstoklet,
     1     mu,istress,mu,nu,ns1,targinfo,nt1,pot,pre,
     2     grad,thresh)
      pot(1) = pot(1)/(2*pi)
      pot(2) = pot(2)/(2*pi)

      do j = 1,2
         do i = 1,2
            ipars(1) = i
            ipars(2) = j
            call st2d_dlp(srcinfo,2,targinfo,0,dpars,0,zpars,
     1           2,ipars,g1)
            if (abs(g1-gmat(i,j)) .gt. 1d-14) ipass=0
         enddo
      enddo
      
      if (sqrt( (pot(1)-u(1))**2 + (pot(2)-u(2))**2) .gt. 1d-14)
     1     ipass=0

c
c     STRAC tests 
c

      call st2d_strac_vec(4,srcinfo,8,targinfo,
     1     0,dpars,0,zpars,0,ipars,tracmat)
      

      trac(1) = tracmat(1,1)*mu(1)+tracmat(1,2)*mu(2)
      trac(2) = tracmat(2,1)*mu(1)+tracmat(2,2)*mu(2)

      thresh = 1d-15
      nd1 = 1
      ns1 = 1
      nt1 = 1
      pot(1) = 0
      pot(2) = 0
      pre = 0
      grad(1,1) = 0
      grad(2,1) = 0
      grad(1,2) = 0
      grad(2,2) = 0
      call st2ddirectstokg(nd1,srcinfo,mu,ns1,targinfo,nt1,
     1     pot,pre,grad,thresh)

      e(1,1) = grad(1,1)
      e(2,1) = (grad(1,2)+grad(2,1))/2
      e(1,2) = (grad(1,2)+grad(2,1))/2
      e(2,2) = grad(2,2)

      nu(1) = targinfo(7)
      nu(2) = targinfo(8)
      trac2(1) = (-pre*nu(1)+2*e(1,1)*nu(1)+2*e(1,2)*nu(2))/(pi*2)
      trac2(2) = (-pre*nu(2)+2*e(2,1)*nu(1)+2*e(2,2)*nu(2))/(pi*2)

      do j = 1,2
         do i = 1,2
            ipars(1) = i
            ipars(2) = j
            call st2d_strac(srcinfo,8,targinfo,0,dpars,0,zpars,
     1           2,ipars,g1)
            if (abs(g1-tracmat(i,j)) .gt. 1d-14) ipass=0
         enddo
      enddo
      
      if (sqrt( (trac(1)-trac2(1))**2
     1     + (trac(2)-trac2(2))**2) .gt. 1d-14) ipass=0

      
c
c     COMB tests
c

      dpars(1) = 0.2d0
      dpars(2) = 0.3d0
      call st2d_comb_vec(4,srcinfo,2,targinfo,2,dpars,0,zpars,
     1     0,ipars,gmat)

      u(1) = gmat(1,1)*mu(1) + gmat(1,2)*mu(2)
      u(2) = gmat(2,1)*mu(1) + gmat(2,2)*mu(2)

      thresh = 1d-15
      nd1 = 1
      ns1 = 1
      nt1 = 1
      nu(1) = srcinfo(7)
      nu(2) = srcinfo(8)
      ifstoklet = 1
      istress=1
      pot(1)=0
      pot(2)=0
      pre = 0
      grad(1,1) = 0
      grad(2,1) = 0
      grad(1,2) = 0
      grad(2,2) = 0
      stoklet(1) = dpars(1)*mu(1)
      stoklet(2) = dpars(1)*mu(2)
      strslet(1) = dpars(2)*mu(1)
      strslet(2) = dpars(2)*mu(2)
      call st2ddirectstokstrsg(nd1,srcinfo,ifstoklet,
     1     stoklet,istress,strslet,nu,ns1,targinfo,nt1,pot,pre,
     2     grad,thresh)
      pot(1) = pot(1)/(2*pi)
      pot(2) = pot(2)/(2*pi)
      
      do j = 1,2
         do i = 1,2
            ipars(1) = i
            ipars(2) = j
            call st2d_comb(srcinfo,2,targinfo,2,dpars,0,zpars,
     1           2,ipars,g1)
            call st2d_slp(srcinfo,2,targinfo,0,dpars,0,zpars,
     1           2,ipars,s1)
            call st2d_dlp(srcinfo,2,targinfo,0,dpars,0,zpars,
     1           2,ipars,d1)
            if (abs(g1-gmat(i,j)) .gt. 1d-14) ipass=0
            if (abs(g1-s1*dpars(1)-d1*dpars(2)) .gt. 1d-14) ipass=0
         enddo
      enddo

      if (sqrt( (pot(1)-u(1))**2 + (pot(2)-u(2))**2) .gt. 1d-14)
     1     ipass=0
      
      return
      end
      
