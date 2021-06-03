      subroutine test_dchunkints(isuc)

      implicit real *8 (a-h,o-z)
      real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8 xytarg(2,100)
      real *8, allocatable :: dintvals(:,:),dintvalsex(:,:)
      integer ipars
      real *8 dpars
      complex *16 zpars

      external dlslp

      done = 1
      pi = atan(done)*4

      call prini(6,13)
      nsuccess = 0

      norder = 8
      allocate(srccoefs(6,norder),srcvals(8,norder))
      allocate(ts(norder),wts(norder),umat(norder,norder),
     1   vmat(norder,norder))  
c
c   get legendre nodes, weights and vals2coefs and coefs2vals 
c   matrices
c
      itype = 2
      call legeexps(itype,norder,ts,umat,vmat,wts)

      xshift = 1.1d0
      yshift = 0.3d0

      do i=1,norder
        srcvals(1,i) = ts(i) + xshift
        srcvals(2,i) = yshift
        srcvals(3,i) = 1
        srcvals(4,i) = 0
        srcvals(5,i) = 0
        srcvals(6,i) = 0
        srcvals(7,i) = 0
        srcvals(8,i) = 1

        do j=1,6
          srccoefs(j,i) = 0
        enddo
      enddo

      do i=1,norder
        do j=1,6
          do l=1,norder
            srccoefs(j,i) = srccoefs(j,i) + umat(i,l)*srcvals(j,l)
          enddo
        enddo
      enddo

      ntarg = 2
      xytarg(1,1) = -1.1d0 + xshift
      xytarg(2,1) = 1.0d0/11.0d0 + yshift
      
      xytarg(1,2) = 0.3d0 + xshift
      xytarg(2,2) = -1.0d0/20.0d0 + yshift

      eps = 1.0d-12
      ndtarg = 2
      ndd = 0
      ndz = 0
      ndi = 0
      nporder = 3
      nqorder = 12
      nchmax = 100000
      allocate(dintvals(nporder,ntarg),dintvalsex(nporder,ntarg))

      call dchunkints_lege(eps,norder,srccoefs,ndtarg,ntarg,xytarg,
     1  nporder,dlslp,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nchmax,
     2  dintvals)

      dintvalsex(1,1) = -0.02812024230649564d0 
      dintvalsex(2,1) =  0.1199608346899685d0
      dintvalsex(3,1) =  -0.02749523721067808d0

      dintvalsex(1,2) =  -0.2791999119798363d0
      dintvalsex(2,2) =  -0.08532857831208408d0
      dintvalsex(3,2) = 0.06982115435807171d0

      call prin2('dintvalsex=*',dintvalsex(1:3,1:2),6)
      call prin2('dintvals=*',dintvals(1:3,1:2),6)
      
      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          ra = ra + dintvalsex(j,i)**2
          erra = erra + (dintvals(j,i)-dintvalsex(j,i))**2
        enddo
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in computing integrals=*',erra,1)
      isuc = 0
      if(abs(erra).lt.eps) isuc = 1

      return
      end




      subroutine dlslp(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(6), y(ndtarg),dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 f

      done = 1
      pi = atan(done)*4
      
      r2 = (x(1) - y(1))**2 + (x(2)-y(2))**2
      f = log(r2)/4/pi

      return
      end
     

