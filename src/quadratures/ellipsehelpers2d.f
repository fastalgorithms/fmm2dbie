

      subroutine ellipse_nearfield2d_getnovers(eps,rho,npolyfac,
     1     nch,norders,ising,novers,ier)
c
c     given bernstein ellipse parameter, precision, and number
c     of polynomials to integrate, returns required order for
c     Gauss Legendre integration rule to accurately integrate 
c     
c     \int_{-1}^1 log|w-t| t^i + t^j \, dt
c
c     where w is outside of the Bernstein ellipse
c     
c     Input:
c
c     eps - real *8, desired precision (smallest: 1d-15)
c     rho - real *8, Bernstein ellipse parameter
c     nch - integer number of chunks
c     norders - integer (nch) array, order of each chunk
c     npolyfac - integer, norders(i)*npolyfac is the
c              number of polynomials to integrate on chunk i
c     ising - integer, type of singularity. currently ignored.
c
c     Output:
c
c     novers - integer (nch), the order of the rule to use
c     ier - error flag
c     ier = 0, normal execution
c

      implicit  real *8 (a-h,o-z)
      integer :: norders(nch), novers(nch)
      integer nlegtab(10000), iprecs(10), npolys(40)
      real *8 :: rhos(10)

      ier=0

      call loadellipseinfo_mem(nrho,nnpoly,nprec)

      call loadellipseinfo(rhos,npolys,iprecs,nlegtab)

      iprec=1
      do i = 1,nprec
         iprec=i
         if (eps .ge. (1d-1)**iprecs(i)) exit
         if (i .eq. nprec) ier = 1
      enddo
      irho=1
      do i = nrho,1,-1
         irho=i
         if (rhos(i) .le. rho) exit
         if (i .eq. 1) ier = 4
      enddo

      do j = 1,nch
         npoly=npolyfac*norders(j)
         ipoly=1
         do i = 1,nnpoly
            ipoly = i
            if (npolys(i) .ge. npoly) exit
            if (i .eq. nnpoly) ier = 2
         enddo
         call ellipse_nearfield2d_readiarr3entry(irho,ipoly,iprec,
     1        nrho,nnpoly,nprec,nlegtab,novers(j))
      enddo
      
      return
      end

      subroutine ellipse_nearfield2d_readiarr3entry(i,j,k,l,m,n,imat,
     1     iout)
      implicit real *8 (a-h,o-z)
      integer :: imat(l,m,n)

      iout = imat(i,j,k)
      
      return
      end


      subroutine ellipse_nearfield2d_definerects(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,rho,rects)
c
c     given a bernstein ellipse parameter, this routine
c     defines a set of rectangles, each of which covers
c     the image of the bernstein ellipse of parameter rho
c     under the chart for a panel in the discretization
c     
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: nch,norders(nch),npts
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts)
      real *8, intent(out) :: rects(2,4,nch)

      complex *16, allocatable :: ellmap(:,:), zcoef(:)
      real *8, allocatable :: ellpts(:,:), cvxpoly(:,:)
      real *8, allocatable :: ellptsort(:,:), ellptwork(:,:)
      complex *16 :: ima, z
      data ima /(0d0,1d0)/

      pi = 4*atan(1d0)
      
      nordermax=0
      do i = 1,nch
         nordermax=max(nordermax,norders(i))
      enddo
      nm1 = nordermax-1
      nell = max(2*nordermax+2,18)

      allocate(ellmap(nordermax,nell),zcoef(nordermax))
      
      do i = 1,nell
         tt = (i-1)*2*pi/nell
         z = rho*exp(ima*tt)
         z = (z+1/z)/2
         call legepolz(z,nm1,ellmap(1,i))
      enddo

      nellp1=nell+1
      allocate(ellpts(2,nell),ellptsort(2,nell),ellptwork(2,nell),
     1     cvxpoly(2,nellp1))

      do i = 1,nch
         norder=norders(i)
         do j = 1,norder
            ii = ixys(i)+j-1
            zcoef(j) = srccoefs(1,ii)+ima*srccoefs(2,ii)
         enddo
         call ellipse_nearfield2d_zmattv(nordermax,norder,nell,
     1        ellmap,zcoef,ellpts)
         call ellipse_nearfield2d_mergesort(nell,ellpts,
     1        ellptsort,ellptwork)

         call convhull2d(nell,ellptsort,nhull,cvxpoly)
         nhullm1=nhull-1
         
         call boundingrect_cvxpolygon(cvxpoly,nhullm1,
     1        rects(1,1,i))
         
      enddo
      
      return
      end


      subroutine ellipse_nearfield2d_mergesort(n,pts,ptsort,ptwork)
c
c     this routine sorts 2d points by x coord (then y if a tie)
c      
      implicit real *8 (a-h,o-z)
      real *8 :: pts(2,n), ptsort(2,n), ptwork(2,n)
      real *8 :: ptl(2), ptr(2)
      maxiter=200
      
      nn = 1

      do j = 1,n
         ptsort(1,j)=pts(1,j)
         ptsort(2,j)=pts(2,j)
      enddo
      
      do i = 1,maxiter
         nn2 = 2*nn
         do j = 1,n,nn2
            ileftend = min(j+nn-1,n)
            irightend = min(j+nn2-1,n)
            ileft = j
            iright = j+nn
            nchunk = irightend-ileft+1
            
            do k = 1,nchunk
               if (ileft .le. ileftend .and.
     1              iright .le. irightend) then
                  ptl(1) = ptsort(1,ileft)
                  ptl(2) = ptsort(2,ileft)
                  ptr(1) = ptsort(1,iright)
                  ptr(2) = ptsort(2,iright)
                  
                  if (ptl(1) .lt. ptr(1) .or.
     1      (ptl(1) .eq. ptr(1) .and. ptl(2) .le. ptr(2)) ) then
                     ptwork(1,j+k-1)=ptl(1)
                     ptwork(2,j+k-1)=ptl(2)
                     ileft = ileft+1
                     if (ileft .le. ileftend) then
                        ptl(1)=ptsort(1,ileft)
                        ptl(2)=ptsort(2,ileft)
                     endif
                  else  
                     ptwork(1,j+k-1)=ptr(1)
                     ptwork(2,j+k-1)=ptr(2)
                     iright = iright+1
                     if (iright .le. irightend) then
                        ptr(1)=ptsort(1,iright)
                        ptr(2)=ptsort(2,iright)
                     endif
                  endif
               else if (ileft .le. ileftend) then
                  ptl(1) = ptsort(1,ileft)
                  ptl(2) = ptsort(2,ileft)
                  
                  ptwork(1,j+k-1)=ptl(1)
                  ptwork(2,j+k-1)=ptl(2)
                  ileft = ileft+1
                  if (ileft .le. ileftend) then
                     ptl(1)=ptsort(1,ileft)
                     ptl(2)=ptsort(2,ileft)
                  endif
               else if (iright .le. irightend) then
                  ptr(1) = ptsort(1,iright)
                  ptr(2) = ptsort(2,iright)
                  
                  ptwork(1,j+k-1)=ptr(1)
                  ptwork(2,j+k-1)=ptr(2)
                  iright = iright+1
                  if (iright .le. irightend) then
                     ptr(1)=ptsort(1,iright)
                     ptr(2)=ptsort(2,iright)
                  endif
               endif               

            enddo

            do k = 1,nchunk
               ptsort(1,j+k-1)=ptwork(1,j+k-1)
               ptsort(2,j+k-1)=ptwork(2,j+k-1)
            enddo
         enddo
         if (nn .gt. n) exit
         nn=nn2
      enddo

      return
      end
      
      
      subroutine ellipse_nearfield2d_zmattv(lda,m,n,a,x,y)
c
c     returns y = a^T * x
c
c     input:
c
c     lda - integer, leading dimension of a array
c     m - integer, number of rows in a, entries in x
c     n - integer, number of cols in a, entries in y
c     a - complex *16 array, lda is leading dimension
c     x - complex *16 array, length m
c
c     output:
c
c     y - complex *16 array, length n
c     

      implicit real *8 (a-h,o-z)
      complex *16 :: a(lda,*), x(*), y(*)

      do i = 1,n
         y(i)=0
         do j = 1,m
            y(i) = y(i)+a(j,i)*x(j)
         enddo
      enddo

      return
      end
