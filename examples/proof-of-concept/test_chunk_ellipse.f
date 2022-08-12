      implicit real *8 (a-h,o-z)
      real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: srcvals(:,:,:),srccoefs(:,:,:),ab(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),adjs(:,:)
      
      complex *16, allocatable :: zcoefs(:,:), zpols(:)
      complex *16, allocatable :: zel(:,:), zcurv(:,:) 
      complex *16 ima, z, zf
      real *8 :: ctr(2), vecs(2,2), ds(2,2), zrect(2,5), dpars(2)
      integer :: ipars(3)
      integer, allocatable :: nlegtab(:,:,:), npolys(:), iprecs(:)
      real *8, allocatable :: rhos(:), xel(:,:), yel(:,:), zgrid(:,:,:)
      real *8, allocatable :: zval(:), poto(:), potov(:), wtso(:)
      real *8, allocatable :: wtsov(:), chargeso(:), chargesov(:)
      real *8, allocatable :: dipveco(:,:), dipvecov(:,:), dipstro(:)
      real *8, allocatable :: dipstrov(:), sigmaover(:) 
      real *8, allocatable :: sourceso(:,:), sourcesov(:,:)
      real *8, allocatable :: srco(:,:), srcov(:,:)
      real *8, allocatable :: sigmaover_very(:), sigma(:,:),to(:)
      real *8, allocatable :: tk(:), tov(:)      
      integer, allocatable :: novers(:), novers_very(:), ixyso(:)
      integer, allocatable :: ixysov(:)
      real *8 :: xrect(10), yrect(10)
      data ima /(0d0,1d0)/
      external fstarn_simple
      
      done = 1
      pi = atan(done)*4

      call prini(6,13)
      nsuccess = 0

      eps = 1d-6
      norder = 16
      npoly = 2*norder
      nover_very = 300
      rho=1.2d0
      

      call prin2('rho *',rho,1)

      call loadellipseinfo_mem(nrho,nnpoly,nprec)
      allocate(nlegtab(nrho,nnpoly,nprec),rhos(nrho),npolys(nnpoly),
     1     iprecs(nprec))

      call loadellipseinfo(rhos,npolys,iprecs,nlegtab)

      do i = nrho,1,-1
         irho=i
         if (rhos(i) .le. rho) exit
      enddo

      do i = 1,nnpoly
         inpoly=i
         if (npolys(i) .ge. npoly) exit
      enddo

      do i = 1,nprec
         iprec = i
         if (10d0**(-iprecs(i)) .le. eps) exit
      enddo

      write(*,*) irho, inpoly, iprec

      nover = nlegtab(irho,inpoly,iprec)

      call prinf('nover *',nover,1)
      
      nchmax = 200000
      k = norder
      npts_max = nchmax*k
      allocate(srcvals(8,k,nchmax),srccoefs(6,k,nchmax),ab(2,nchmax))
      allocate(norders(nchmax),ixys(nchmax+1),iptype(nchmax))
      allocate(adjs(2,nchmax)) 
      allocate(zcoefs(k,nchmax),zcurv(k,nchmax),zpols(k))     

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

! parameters sent to routine defining the curve (fstarn_simple)
      dpars(1) = 1.0d0
      dpars(2) = 0.25d0
      ipars(1) = 5

! discretization paratmers
      novergeo = 0
      nch = 0
      ier = 0
      epsc = eps*100
      epsc = 1d-5
      call chunkfunc_guru(epsc,rlmax,ifclosed,irefinel,irefiner, 
     1     rlmaxe,ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi, 
     2     ipars,novergeo,k,nchmax,nch,norders,ixys,iptype,npts, 
     3     srcvals,srccoefs,ab,adjs,ier)

      call prinf('nch *',nch,1)
      call prinf('ier *',ier,1)
      if (ier .ne. 0) stop

      
      do j = 1,nch
         do i = 1,k
            zcoefs(i,j) = srccoefs(1,i,j) + ima*srccoefs(2,i,j)
            zcurv(i,j) = srcvals(1,i,j) + ima*srcvals(2,i,j)
         enddo
      enddo

      call prin2('zcoefs *',zcoefs,2*k+2)

      jimage = 9


      nm1 = norder-1
      nplot = max(2*norder+2,16)

      allocate(zel(nplot,nch))

      do jj = 1,nch
         do i = 1,nplot
            tt = (i-1)*2*pi/nplot
            z = rho*exp(ima*tt)
            z = (z+1/z)/2
            call legepolz(z,nm1,zpols)
            zel(i,jj)=0
            do j = 1,k
               zel(i,jj)= zel(i,jj) + zpols(j)*zcoefs(j,jj)
            enddo
         enddo
      enddo

      nplot_tot = nplot*nch

      jrect= jimage
      call get_rect(srcvals(1,1,jrect),norder,zel(1,jrect),
     1     nplot,ctr,vecs,ds)

      zrect(1,1) = ctr(1)+vecs(1,1)*ds(1,1)+vecs(1,2)*ds(1,2)
      zrect(2,1) = ctr(2)+vecs(2,1)*ds(1,1)+vecs(2,2)*ds(1,2)
      zrect(1,2) = ctr(1)+vecs(1,1)*ds(1,1)+vecs(1,2)*ds(2,2)
      zrect(2,2) = ctr(2)+vecs(2,1)*ds(1,1)+vecs(2,2)*ds(2,2)
      zrect(1,3) = ctr(1)+vecs(1,1)*ds(2,1)+vecs(1,2)*ds(2,2)
      zrect(2,3) = ctr(2)+vecs(2,1)*ds(2,1)+vecs(2,2)*ds(2,2)
      zrect(1,4) = ctr(1)+vecs(1,1)*ds(2,1)+vecs(1,2)*ds(1,2)
      zrect(2,4) = ctr(2)+vecs(2,1)*ds(2,1)+vecs(2,2)*ds(1,2)
      zrect(1,5) = zrect(1,1)
      zrect(2,5) = zrect(2,1)

      do i = 1,5
         xrect(i)=zrect(1,i)
         yrect(i)=zrect(2,i)
      enddo

      x0 = 1d300
      x1 = -1d300
      y0=x0
      y1=x1
      do i = 1,k
         x0 = min(x0,srcvals(1,i,jimage))
         x1 = max(x1,srcvals(1,i,jimage))
         y0 = min(y0,srcvals(2,i,jimage))
         y1 = max(y1,srcvals(2,i,jimage))
      enddo
      dx = x1-x0
      dy = y1-y0
      dd = max(dx,dy)
      xmid = (x0+x1)/2
      ymid = (y0+y1)/2
      fac = 0.75d0
      x0 = xmid-dd*fac
      y0 = ymid-dd*fac
      x1 = xmid+dd*fac
      y1 = ymid+dd*fac

c     set up a grid

      nx = 400
      ny = nx
c      y0 = -1.5
c      y1 = 1.5
c      x0 = -1.5
c      x1 = 1.5
      xl = (x1-x0)/nx
      yl = (y1-y0)/ny

      allocate(zgrid(2,ny,nx))

c     
      
      do j = 1,nx
         x = (j-0.5d0)*xl + x0
         do i = 1,ny
            y = -(i-0.5d0)*yl + y1

            zgrid(1,i,j) = x
            zgrid(2,i,j) = y
         enddo
      enddo

      write(*,*) 'x0,..', x0,x1,y0,y1

      allocate(xel(nplot,nch),yel(nplot,nch))
      ii=0
      do jj = 1,nch
         do i = 1,nplot
            ii = ii +1
            xel(i,jj) = real(zel(i,jj))
            yel(i,jj) = -real(ima*zel(i,jj))
         enddo
      enddo

c     compare oversampled to very oversampled rule...

c     set up some charge density  

      allocate(sigma(k,nch))
      
      do j = 1,nch
         do i = 1,k
            sigma(i,j) = cos(300*srcvals(1,i,j)-30*srcvals(2,i,j))
         enddo
      enddo

c     oversample
      
      allocate(novers(nch),novers_very(nch),ixyso(nch+1),ixysov(nch+1))
      ixyso(1)=1
      ixysov(1)=1
      do i = 1,nch
         novers(i)=nover
         novers_very(i)=nover_very
         ixyso(i+1)=ixyso(i)+nover
         ixysov(i+1)=ixysov(i)+nover_very
      enddo

      nso = nch*nover
      nsov = nch*nover_very

      allocate(srco(8,nso),srcov(8,nsov))

      write(*,*) 'oversampling geometry ..'
      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,novers,ixyso,nso,srco)
      write(*,*) 'oversampling geometry (a lot)..'      
      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,novers_very,ixysov,nsov,srcov)

      allocate(wtso(nso),wtsov(nsov))

      write(*,*) 'getting oversampled weights'
      call get_qwts2d(nch,novers,ixyso,iptype,nso,
     1     srco,wtso)
      write(*,*) 'getting very oversampled weights'      
      call get_qwts2d(nch,novers_very,ixysov,iptype,nsov,
     1     srcov,wtsov)

      
      allocate(sigmaover(nso),sigmaover_very(nsov))
      call oversample_fun_curv2d(1,nch,norders,ixys,iptype, 
     1     npts,sigma,novers,ixyso,nso,sigmaover)
      call oversample_fun_curv2d(1,nch,norders,ixys,iptype, 
     1     npts,sigma,novers_very,ixysov,nsov,sigmaover_very)

c     call fmm
      
      allocate(sourceso(2,nso),sourcesov(2,nsov))
      allocate(chargeso(nso),chargesov(nsov)) 
      allocate(dipstro(nso),dipstrov(nsov)) 
      allocate(dipveco(2,nso),dipvecov(2,nsov)) 
      
      alpha = 1
      beta = 0
      pi2 = pi*2
      
      do i=1,nso
         sourceso(1,i) = srco(1,i)
         sourceso(2,i) = srco(2,i)

         chargeso(i) = -sigmaover(i)*wtso(i)*alpha/pi2
         dipstro(i) = -sigmaover(i)*wtso(i)*beta/pi2
         dipveco(1,i) = srco(7,i)
         dipveco(2,i) = srco(8,i)
      enddo
      
      do i=1,nsov
         sourcesov(1,i) = srcov(1,i)
         sourcesov(2,i) = srcov(2,i)

         chargesov(i) = -sigmaover_very(i)*wtsov(i)*alpha/pi2
         dipstrov(i) = -sigmaover_very(i)*wtsov(i)*beta/pi2
         dipvecov(1,i) = srcov(7,i)
         dipvecov(2,i) = srcov(8,i)
      enddo

      ifcharge=0
      ifdipole=0
      if (abs(alpha) .gt. 1d-15) ifcharge = 1
      if (abs(beta) .gt. 1d-15) ifbeta = 1
      
      epsfmm=max(eps/100,1d-15)
      epsfmm=1d-15
      ifpgh = 0
      ifpghtarg = 1
      ntarg = nx*ny
      allocate(poto(ntarg),potov(ntarg))

      call prinf('ifcharge *',ifcharge,1)
      call prin2('charges *',chargeso,100)      
      
      write(*,*) 'calling oversampled fmm ', nso, ntarg

      nd=1
      call rfmm2d(nd,epsfmm,nso,sourceso,ifcharge,chargeso,
     1     ifdipole,dipstro,dipveco,iper,ifpgh,tmp,tmp,tmp,ntarg,
     1     zgrid,ifpghtarg,poto,tmp,tmp,ier)

      write(*,*) 'calling very oversampled fmm ', nsov, ntarg
      
      call rfmm2d(nd,epsfmm,nsov,sourcesov,ifcharge,chargesov,
     1     ifdipole,dipstrov,dipvecov,iper,ifpgh,tmp,tmp,tmp,ntarg,
     1     zgrid,ifpghtarg,potov,tmp,tmp,ier)


      call prin2('potov *',potov,20)
      call prin2('poto *',poto,20)
      

      allocate(zval(ntarg))
      
      do i = 1,ntarg
         derr = abs(poto(i)-potov(i))/max(abs(potov(i)),1d0)
         zval(i) = log(derr)/log(1d1)
      enddo
         

      
      iw = 1
      icmp=6
c      call pyimage6(iw,ny,nx,zval,xel,yel,nplot_tot,x0,x1,
c     1     y0,y1,icmp)


      write(*,*) 'x0,..', x0,x1,y0,y1
      nrect=5
      call pyimage7(iw,ny,nx,zval,xel(1,jimage),yel(1,jimage),nplot,
     1     xrect,yrect,nrect,x0,x1, y0,y1,icmp)

      allocate(to(nover),tov(nover_very),tk(k))
      itype=0
      call legeexps(itype,k,tk,tmp,tmp,tmp)
      call legeexps(itype,nover,to,tmp,tmp,tmp)
      call legeexps(itype,nover_very,tov,tmp,tmp,tmp)

c      jj = 15
c      call pyplot3(iw,to,sigmaover(ixyso(jj)),nover,1,tov,
c     1     sigmaover_very(ixysov(jj)),
c     1     nover_very,3,tk,sigma(1,jj),k,1,'compare *')
      
c      
c      iw=1
c      call zpyplot3(iw,zel,nplot_tot,1,zcurv,npts,3,zrect,5,3,
c     1     'curve and ellipse *')
c
      stop
      end

      subroutine get_circ(srcvals,norder,relpts,nrel,ctr,rad)
      implicit real *8 (a-h,o-z)
      dimension srcvals(8,norder), relpts(2,nrel),ctr(2)

      ctr(1) = 0
      ctr(2) = 0
      do i = 1,norder
         ctr(1) = ctr(1)+srcvals(1,i)
         ctr(2) = ctr(2)+srcvals(2,i)
      enddo

      ctr(1) = ctr(1)/norder
      ctr(2) = ctr(2)/norder
      
      rmax = 0
      do i = 1,nrel
         x = relpts(1,i)-ctr(1)
         y = relpts(2,i)-ctr(2)
         rmax = max(rmax,x**2+y**2)
      enddo

      rad = sqrt(rmax)

      return
      end
      
      subroutine get_rect(srcvals,norder,relpts,nrel,ctr,vecs,ds)
      implicit real *8 (a-h,o-z)
      dimension srcvals(8,norder), relpts(2,nrel),ctr(2),vecs(2,2)
      dimension ds(2,2)

      real *8 :: dnave(2), dtave(2)

      dnave(1) = 0
      dnave(2) = 0
      ctr(1) = 0
      ctr(2) = 0
      do i = 1,norder
         dnave(1) = dnave(1)+srcvals(7,i)
         dnave(2) = dnave(2)+srcvals(8,i)
         ctr(1) = ctr(1)+srcvals(1,i)
         ctr(2) = ctr(2)+srcvals(2,i)
      enddo

      ctr(1) = ctr(1)/norder
      ctr(2) = ctr(2)/norder
      
      dd = sqrt(dnave(1)**2 + dnave(2)**2)
      dnave(1) = dnave(1)/dd
      dnave(2) = dnave(2)/dd

      dtave(1) = -dnave(2)
      dtave(2) = dnave(1)

      dnmax = 0
      dnmin = 0
      dtmax = 0
      dtmin = 0

      do i = 1,nrel
         x = relpts(1,i)-ctr(1)
         y = relpts(2,i)-ctr(2)
         dnproj = dnave(1)*x + dnave(2)*y
         dtproj = dtave(1)*x + dtave(2)*y
         dnmax = max(dnmax,dnproj)
         dnmin = min(dnmin,dnproj)
         dtmax = max(dtmax,dtproj)
         dtmin = min(dtmin,dtproj)
      enddo

      vecs(1,1) = dnave(1)
      vecs(2,1) = dnave(2)
      vecs(1,2) = dtave(1)
      vecs(2,2) = dtave(2)

      ds(1,1) = dnmax
      ds(2,1) = dnmin
      ds(1,2) = dtmax
      ds(2,2) = dtmin

      return
      end
      
      subroutine curve_fun(t,srcvals)
      implicit real *8 (a-h,o-z)
      dimension srcvals(8)

      srcvals(1) = t + 2d0
      srcvals(2) = sin(t) + 1.5d0
      srcvals(3) = 1
      srcvals(4) = cos(t)
      srcvals(5) = 0
      srcvals(6) = -sin(t)
      srcvals(7) = srcvals(4)
      srcvals(8) = -srcvals(3)

      ss = sqrt(srcvals(7)**2+srcvals(8)**2)
      srcvals(7) = srcvals(7)/ss
      srcvals(8) = srcvals(8)/ss

      return
      end
      

      subroutine fstarn_simple(t,ndd,dpars,ndz,zpars,ndi,ipars,x, 
     1     y,dxdt,dydt,dxdt2,dydt2)
c     
c     curve parametrization for simple star shaped domain
c     
c     r(t) = r1 + r2*cos(n*t)
c     x(t) = r(t) cos(t)
c     y(t) = r(t) sin(t)
c     
c     where r1 = dpars(1), r2 = dpars(2) and n = ipars(1)
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
     1     ima*(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t)
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
      
