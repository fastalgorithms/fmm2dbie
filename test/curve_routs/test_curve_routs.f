      implicit real *8 (a-h,o-z)

      call prini(6,13)
      
      call test_star3(isuccess)

      return
      end



      subroutine test_star3(isuccess)
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),ab(:,:)
      integer, allocatable :: norders(:),ixys(:),iptype(:),adjs(:,:)
      real *8, allocatable :: xs(:),ys(:)
      real *8, allocatable :: xs2(:),ys2(:)
      real *8 dpars(2)
      integer ipars(1)
      complex *16 zpars

      external fstarn_simple
      

      done = 1
      pi = atan(done)*4


      nchmax = 20000
      k = 12
      npts_max = nchmax*k
      allocate(srcvals(8,npts_max),srccoefs(6,npts_max),ab(2,nchmax))
      allocate(norders(nchmax),ixys(nchmax+1),iptype(nchmax))
      allocate(adjs(2,nchmax))

      irefinel = 0
      irefiner = 0
      rlmax = 1.0d10
      ta = 0
      tb = 2*pi
      ifclosed = 1
      rlmaxe = 0

      ndd = 2
      ndi = 1
      ndz = 0
      
      dpars(1) = 1.0d0
      dpars(2) = 0.5d0

      ipars(1) = 60
      nover = 1
      nch = 0
      ier = 0
      eps = 1.0d-5
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
      allocate(xs(npts),ys(npts))
      do i=1,npts
         xs(i) = srcvals(1,i)
         ys(i) = srcvals(2,i)
      enddo

      allocate(xs2(nch),ys2(nch))

      do i=1,nch
        call fstarn_simple(ab(1,i),ndd,dpars,ndz,zpars,ndi,ipars,
     1     xs2(i),ys2(i),dx,dy,dx2,dy2)

      enddo


      call pyplot2(18,xs2,ys2,nch,1,xs,ys,npts,3,'a*')
      

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
      complex *16 zpars

      rsc = 1.0d0

      n = ipars
      r1 = dpars(1)*rsc
      r2 = dpars(2)*rsc

      r = r1 + r2*cos(n*t)
      drdt = -r2*n*sin(n*t)
      d2rdt2 = -r2*n*n*cos(n*t)

      x = r*cos(t)
      y = r*sin(t)
      
      dxdt = drdt*cos(t) - r*sin(t)
      dydt = drdt*sin(t) + r*cos(t)
      
      dxdt2 = d2rdt2*cos(t) - 2*drdt*sin(t) - r*cos(t) 
      dydt2 = d2rdt2*sin(t) + 2*drdt*cos(t) - r*sin(t) 
      
      return
      end



