      implicit real *8 (a-h,o-z)

      call prini(6,13)
      
      call test_chunknear1(isuccess1)
      call test_findnearchunktarg_id_ts(isuccess2)

      call prinf('isuccess1 *',isuccess1,1)
      call prinf('isuccess2 *',isuccess2,1)      

      return
      end



      subroutine test_chunknear1(isuccess)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),ab(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),adjs(:,:)
      real *8, allocatable :: xs(:),ys(:)
      real *8, allocatable :: xs2(:),ys2(:),targs(:,:),ttest(:)
      real *8, allocatable :: xtest(:),ytest(:)      
      real *8, allocatable :: tinits(:),pjcoefs1(:),pjcoefs2(:)
      real *8, allocatable :: xcoefs(:),ycoefs(:),pexp0(:),ts_targ(:)
      real *8, allocatable :: dist2_targ(:),cms(:,:),rads(:)     
      real *8 dpars(2), tt(30), ww(30), uu, vv
      integer ipars(1), ichtest(100)
      integer, allocatable :: idx(:)
      complex *16 zpars

      external fstarn_simple
      

      done = 1
      pi = atan(done)*4


      nchmax = 200000
      k = 14
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
      
      dpars(1) = 1.0d0
      dpars(2) = 0.5d0

      ipars(1) = 5
      nover = 1
      nch = 0
      ier = 0
      eps = 1.0d-2
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1  ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi,ipars,nover,
     2  k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,adjs,
     3  ier)
      
      call prinf('ier=*',ier,1)
      if(ier.ne.0) then
        call prinf('failed to execute*',ier,1)
        stop
      endif
      call prinf('nch=*',nch,1)

      allocate(cms(2,nch),rads(nch))
      call get_centroid_rads2d(nch,norders,ixys,iptype,npts,srccoefs,
     1     srcvals,cms,rads)

      

      iseed = 12345
      d = hkrand(iseed)
      nchtest = 10
      do i = 1,nchtest
         ichtest(i) = hkrand(0)*nch
         ichtest(i) = max(ichtest(i),1)
         ichtest(i) = min(ichtest(i),nch)
      enddo

      call prinf('ichtest *',ichtest,nchtest)


      ngrid = 100
      
      nt = ngrid*ngrid

      allocate(targs(2,nt))
      allocate(ts_targ(nt),dist2_targ(nt))
      allocate(idx(nt))

      ntinit = 2*k+2
      ncoef = k+3
      allocate(tinits(ntinit),pjcoefs1(ncoef),pjcoefs2(ncoef),
     1     xcoefs(k),ycoefs(k),pexp0(k))

      ntest = 10000
      allocate(ttest(ntest),xtest(ntest),ytest(ntest))

      do i = 1,k
         xcoefs(i) = srccoefs(1,i)
         ycoefs(i) = srccoefs(2,i)
      enddo

      call prin2('xcoefs *',xcoefs,k)
      call prin2('ycoefs *',ycoefs,k)

      
      do iii = 1,nchtest
         ich = ichtest(iii)
         
c     test the near chunk routine on the first chunk
         
         xmin = srcvals(1,1+(ich-1)*k)
         xmax = srcvals(1,1+(ich-1)*k)
         ymin = srcvals(2,1+(ich-1)*k)
         ymax = srcvals(2,1+(ich-1)*k)
         
         do i = 1,k
            x = srcvals(1,i+(ich-1)*k)
            y = srcvals(2,i+(ich-1)*k)
            if (x .lt. xmin) xmin = x
            if (x .gt. xmax) xmax = x
            if (y .lt. ymin) ymin = y
            if (y .gt. ymax) ymax = y
         enddo

         dx = xmax-xmin
         xmin = xmin - 4*dx
         xmax = xmax + 4*dx
         dy = ymax-ymin
         ymin = ymin - 4*dy
         ymax = ymax + 4*dy


         dx = (xmax-xmin)/ngrid
         dy = (ymax-ymin)/ngrid
         
         ii = 0
         do i = 1,ngrid
            do j = 1,ngrid
               ii = ii + 1
               targs(1,ii) = xmin+dx/2+(i-1)*dx
               targs(2,ii) = ymin+dy/2+(j-1)*dy
            enddo
         enddo

         do i = 1,nt
            idx(i) = i
         enddo


         itype=0
         ntinit2 = ntinit-2
         call legeexps(itype,ntinit2,tinits(2),u,v,w)
         tinits(1) = -1
         tinits(ntinit) = 1

         ndeg = k-1
         do i=1,k
            pexp0(i) = 0
         enddo
         x0 = 0
         call legeexe2(x0,val0,pexp0,ndeg,pjcoefs1,pjcoefs2,ncoef)

         do i = 1,k
            xcoefs(i) = srccoefs(1,i+(ich-1)*k)
            ycoefs(i) = srccoefs(2,i+(ich-1)*k)
         enddo

         chlen = rads(ich)
         
         ndt = 2
         norderloc = k
         ifprint = 1
         call cpu_time(t1)
         call chunknear1(norderloc,xcoefs,ycoefs,chlen,ntinit,tinits,
     1        ncoef,pjcoefs1,pjcoefs2,ndt,nt,targs,nt,idx,ts_targ,
     2        dist2_targ,ifprint)
         call cpu_time(t2)
         write(*,*) iii, nt/(t2-t1)

c     test against shifted problem (that's what the routine
c     optimizes for)

         do i = 1,k
            xcoefs(i) = srccoefs(1,i+(ich-1)*k)
            ycoefs(i) = srccoefs(2,i+(ich-1)*k)
         enddo

         xc1 = xcoefs(1)
         yc1 = ycoefs(1)
c         xcoefs(1) = 0
c         ycoefs(1) = 0
         

         do i = 1,ntest
            ttest(i) = -1 + (i-1)*2.0d0/(ntest-1)
            call legeexev(ttest(i),xtest(i),xcoefs,ndeg)
            call legeexev(ttest(i),ytest(i),ycoefs,ndeg)
         enddo

         nfail = 0
         do i = 1,nt

c            xt = targs(1,i) - xc1
c            yt = targs(2,i) - yc1
            xt = targs(1,i)
            yt = targs(2,i)
            
            call legeexev(ts_targ(i),xfound,xcoefs,ndeg)
            call legeexev(ts_targ(i),yfound,ycoefs,ndeg)

            df = (xfound-xt)**2 + (yfound-yt)**2

            dtestmin = (xtest(1)-xt)**2 + (ytest(1)-yt)**2
            ibest = 1
            do j = 2,ntest
               dtemp = (xtest(j)-xt)**2 + (ytest(j)-yt)**2
               if (dtemp .lt. dtestmin) then
                  dtestmin = dtemp
                  ibest = j
               endif
            enddo

            if (df-dtestmin .gt. 1d-15*(dtestmin+abs(xc1)+abs(yc1)))
     1           then
               nfail = nfail + 1
               write(*,*) i
               write(*,*) df, dtestmin
               write(*,*) ts_targ(i), ttest(ibest)
               write(*,*) abs(df-dtestmin)/dtestmin,
     1              abs(ts_targ(i) - ttest(ibest))
            endif
            
         enddo

         isuccess = 1
         if (nfail .gt. 0) isuccess = 0

      enddo
      return
      end
c
c
      subroutine test_findnearchunktarg_id_ts(isuccess)
      implicit real *8 (a-h,o-z)
      real *8 :: timeinfo(8)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),ab(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),adjs(:,:)
      real *8, allocatable :: srcrad(:), srcrad1(:), xys(:,:)
      real *8, allocatable :: targs(:,:), ts_targ(:), dist_targ(:)
      real *8, allocatable :: dist_targ2(:), ts_targ2(:)
      integer, allocatable :: ich_id(:), ich_id2(:)
      real *8 dpars(2), tt(30), ww(30), uu, vv
      integer ipars(1), ichtest(100)
      complex *16 zpars

      external fstarn_simple
      

      done = 1
      pi = atan(done)*4


      nchmax = 200000
      k = 14
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
      
      dpars(1) = 1.0d0
      dpars(2) = 0.5d0

      ipars(1) = 5
      nover = 1
      nch = 0
      ier = 0
      eps = 1.0d-13
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1  ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi,ipars,nover,
     2  k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,adjs,
     3  ier)
      
      call prinf('ier=*',ier,1)
      if(ier.ne.0) then
        call prinf('failed to execute*',ier,1)
        stop
      endif
      call prinf('nch=*',nch,1)

      
c     set up a grid 
      
      xmin = srcvals(1,1)
      xmax = srcvals(1,1)
      ymin = srcvals(2,1)
      ymax = srcvals(2,1)
      
      do i = 1,npts
         x = srcvals(1,i)
         y = srcvals(2,i)
         if (x .lt. xmin) xmin = x
         if (x .gt. xmax) xmax = x
         if (y .lt. ymin) ymin = y
         if (y .gt. ymax) ymax = y
      enddo

      dx = xmax-xmin
      xmin = xmin - 4*dx
      xmax = xmax + 4*dx
      dy = ymax-ymin
      ymin = ymin - 4*dy
      ymax = ymax + 4*dy

      ngrid = 100
      
      dx = (xmax-xmin)/ngrid
      dy = (ymax-ymin)/ngrid

      nt = ngrid*ngrid
      ndt = 2
      allocate(targs(ndt,nt))
      
      ii = 0
      do i = 1,ngrid
         do j = 1,ngrid
            ii = ii + 1
            targs(1,ii) = xmin+dx/2+(i-1)*dx
            targs(2,ii) = ymin+dy/2+(j-1)*dy
         enddo
      enddo


      allocate(ich_id(nt),ts_targ(nt),dist_targ(nt),srcrad(npts),
     1     srcrad1(nch),xys(2,nch),ich_id2(nt),dist_targ2(nt),
     2     ts_targ2(nt))

      call get_centroid_rads2d(nch,norders,ixys,iptype,npts,
     1     srccoefs,srcvals,xys,srcrad1)
      
      do i = 1,nch
         istart = ixys(i)
         iend = ixys(i+1)-1
         do j = istart,iend
            srcrad(j) = srcrad1(i)*3.1d0
         enddo
      enddo
      
      norderloc = k
      call cpu_time(t1)
c$    t1=omp_get_wtime()
      call findnearchunktarg_id_ts(nch,norders,ixys,iptype,npts,
     1     srccoefs,srcvals,srcrad,ndt,nt,targs,ich_id,ts_targ,
     2     dist_targ,timeinfo,ier)
      call cpu_time(t2)
c$    t2=omp_get_wtime()

      nnz = timeinfo(3)
      
      spdtree = (nt+npts)/(timeinfo(1))
      spdopt = (nnz)/(timeinfo(2))

      call prinf('after findnearchunktarg_id_ts, ier *',ier,1)
      call prin2('after findnearchunktarg_id_ts, timeinfo *',timeinfo,3)
      call prin2('points per second, findnear basic *',spdtree,1)
      call prin2('points per second, optimization *',spdopt,1)      

      call findnearchunktarg_id_ts_brute(nch,norders,ixys,iptype,
     1     npts,srccoefs,srcvals,ndt,nt,targs,ich_id2,ts_targ2,
     2     dist_targ2)

      ifail = 0
      do i = 1,nt
         if (ich_id(i) .gt. 0 .and. ich_id2(i) .ne. ich_id(i)) then
            if (srcrad(ixys(ich_id2(i))) .gt. dist_targ2(i)) then
               write(*,*)
               write(*,*) 'i ', i
               write(*,*) ich_id(i), ich_id2(i)
               write(*,*) dist_targ(i), dist_targ2(i)
               write(*,*) srcrad(ixys(ich_id(i))),
     1              srcrad(ixys(ich_id2(i)))
               ifail = ifail+1
            endif
         endif
      enddo

      isuccess = 1
      if (ifail .gt. 0) isuccess = 0

      
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

