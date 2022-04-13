      subroutine test_zchunkints(isuc)

      implicit real *8 (a-h,o-z)
      real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:)
      real *8 xytarg(2,100)
      complex *16, allocatable :: zintvals(:,:),zintvalsex(:,:)
      integer ipars
      real *8 dpars
      complex *16 zpars,ima
      data ima/(0.0d0,1.0d0)/

      external zldlp

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

      eps = 1.0d-7
      ndtarg = 2
      ndd = 0
      ndz = 0
      ndi = 0
      nporder = 3
      nqorder = 26
      nchmax = 100000
      allocate(zintvals(nporder,ntarg),zintvalsex(nporder,ntarg))

      call zchunkints_lege(eps,norder,srccoefs,ndtarg,ntarg,xytarg,
     1  nporder,zldlp,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nchmax,
     2  zintvals)

      zintvalsex(1,1) = -0.4367645330548260d0-0.1105413886042094d0*ima  
      zintvalsex(2,1) = 0.1721803173223552d0+0.0818896608232826d0*ima
      zintvalsex(3,1) = -0.07688202898492064d0-0.05636811187599036d0*ima

      zintvalsex(1,2) = 0.0982358058709252d0+0.4825327362729426d0*ima
      zintvalsex(2,2) = -0.2647125076088660d0+0.1398480305883365d0*ima
      zintvalsex(3,2) = -0.1577499290653270d0-0.1584813163010549d0*ima

      call prin2('zintvalsex=*',zintvalsex(1:3,1:2),12)
      call prin2('zintvals=*',zintvals(1:3,1:2),12)
      
      erra = 0
      ra = 0
      do i=1,ntarg
        do j=1,3
          ra = ra + abs(zintvalsex(j,i))**2
          erra = erra + abs(zintvals(j,i)-zintvalsex(j,i))**2
        enddo
      enddo
      erra = sqrt(erra/ra)
      call prin2('error in computing integrals=*',erra,1)

      isuc = 0
      if(abs(erra).lt.eps) isuc = 1

      return
      end




      subroutine zldlp(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
      implicit real *8 (a-h,o-z)
      real *8 x(6), y(ndtarg),dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      complex *16 f,z

      done = 1
      pi = atan(done)*4

      z = dcmplx(y(1)-x(1),y(2)-x(2))
      f = 1.0d0/2/pi/z

      return
      end
     

