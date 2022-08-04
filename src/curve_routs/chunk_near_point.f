c
c
c
c     end of the debugging code and start of the near point code
c     proper
c
c      
c     This file contains 2 user callable routines for finding
c     the nearest chunks for a collection of points
c
c     findnearchunktarg_id_ts - find the nearest chunk (by id)
c         and the parameter t in the chunk discretization 
c         which is closest for each target. The search is 
c         accelerated by only considering chunks for which the
c         target is contained in a user-specified radius 
c         of the boundary points
c     
c     findnearchunktarg_id_ts_brute - find the nearest chunk
c         (by id) and the parameter t in the chunk discretization 
c         which is closest for each target. There is no
c         acceleration so it is O(nt*nch) with a rather large
c         constant. Useful for debugging.
c
c     In addition to these routines you will find the workhorse
c     routines:
c
c     chunknear1 - find the parameter t in the chunk discretization
c          which minimizes the distance between the chunk point
c          and target for each target. Requires precomputed arrays
c     
c     legefdd2 - similar to legefde2 but also computes second
c          derivatives


      subroutine findnearchunktarg_id_ts(nch,norders,ixys,iptype,npts,
     1     srccoefs,srcvals,srcrad,ndt,nt,targets,ich_id,ts_targ,
     2     dist_targ,timeinfo,ier)
c
c
c     Find nearest chunk and point on chunk for each target.
c      
c     This is the faster version in which the search is limited
c     to neighborhoods of each source point.
c
c     Input:
c
c     nch - integer, number of chunks
c     norders - integer(nch), discretization order of each chunk
c     ixys - integer(nch+1), index into srccoefs, srcvals arrays
c     ixys - integer(nch+1)
c       starting location in srccoefs,srcvals array
c       where info for each chunk begins
c     iptype - integer(nch)
c       type of chunk
c       iptype = 1, chunk discretized with Gauss-Legendre
c       nodes and not at a corner
c     npts - integer
c       total number of discretization points
c     srccoefs - real *8 (6,npts)
c       Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c       if(iptype.eq.1) then basis = Legendre polynomials
c     srcvals - real *8 (8,npts)
c     x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization
c       points
c     srcrad - real *8 (npts), suggested search radius for each
c       point
c     ndt - integer, first dimension of targets array
c     nt - integer, second dimension of targets array, number of targs
c     targets - real *8 (ndt,nt), array in which first two
c         entries of targets(:,i) are the coordinates of target i
c     
c     Output:
c
c     NOTE there is nothing which prevents this routine from 
c     identifying a chunk as the closest chunk to the target
c     even though the given target is in fact outside
c     of the radii specified in srcrad (in some interpolated
c     sense). I.e. sufficient closeness is not enforced.
c
c     ich_id - integer(nt), ich_id(i) is the chunk id for chunk
c         nearest to the ith target.
c         if no chunks found using the respective search radii
c         ich_id(i) = -1
c     ts_targ - real *8(nt), ts_targ(i) is the value in [-1,1]
c         of the underlying chunk parameterization of chunk ich_id(i)
c         which is closest to the target. If no close chunk was found
c         ts_targ(i) = 0
c     dist_targ - real *8(nt), dist_targ(i) is the distance
c         from the ith target to the point corresponding to 
c         ts_targ(i) on chunk ich_id(i). If no close chunk found,
c         dist_targ(i) = -1    
c     timeinfo - real *8(3), cpu time (wall time for openmp) of
c         the major tasks 
c                 timeinfo(1) - findnear2d call
c                 timeinfo(2) - time for nnz ts calcs (optimizations)
c                 timeinfo(3) - nnz + 0.1d0
c     ier - integer, error flag. ier .eq. 0 means normal operation
c                 ier .eq. 4, insufficient memory for a smaller
c                          allocation step. check inputs?
c                 ier .lt. 0, insufficient memory, consider
c                          decreasing chunk radii. -ier is the
c                          total number of chunk-target pairs 
c                          for which a distance is to be computed
      implicit none
      integer, intent(in) :: nch,norders(nch),ixys(nch+1),iptype(nch)
      integer, intent(in) :: npts, nt
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),
     1     srcrad(npts), targets(ndt,nt)
      integer, intent(out) :: ich_id(nt), ier
      real *8, intent(out) :: ts_targ(nt), dist_targ(nt), timeinfo(3)

      integer :: ndt, nnz
      real *8, allocatable :: xys(:,:), rads1(:), rads(:), dist22(:)
      real *8, allocatable :: ts_targ2(:), qwts(:), chlens(:)
      
      integer, allocatable :: row_ptr(:), col_ind(:), col_ptr(:)
      integer, allocatable :: row_ind(:), iper(:), iperinv(:)

      real *8 :: pjcoefs1(100), pjcoefs2(100), pexp0(100), xcoefs(100)
      real *8 :: ycoefs(100), tinits(200)

      real *8 :: val0, w(100), u, v, tbest, x0
      real *8 :: t, radmax, d2max, d2, d2best

      real *8 :: t1, t2, chlen
c$    real *8 :: omp_get_wtime

      integer :: ntinit, nt1, norderloc, ndeg, ncoef, k, itype, ii, j
      integer :: istart, iend, ichbest, i, ier1, jj, ifprint
      
      ier = 0

      ier1 = 0
      
c     find chunks sufficiently close to each target
      
      allocate(xys(2,nch),rads1(nch),rads(nch),stat=ier1)
      if (ier1 .ne. 0) then
         ier = 4
         return
      endif

      call get_centroid_rads2d(nch,norders,ixys,iptype,npts,
     1     srccoefs,srcvals,xys,rads1)

      ii = 0
      do i = 1,nch
         radmax = 0
         do j = 1,norders(i)
            ii = ii+1
            radmax = max(radmax,srcrad(ii))
         enddo
         rads(i) = rads1(i) + radmax*1.2d0
      enddo

      call cpu_time(t1)
c$    t1 = omp_get_wtime()       
      call findnear2dmem(xys,nch,rads,ndt,targets,nt,nnz)
      
      allocate(row_ptr(nt+1),col_ind(nnz),stat=ier1)
      if (ier1 .ne. 0) then
         ier = -nnz
         return
      endif
      

      call findnear2d(xys,nch,rads,ndt,targets,nt,row_ptr,
     1     col_ind)

      call cpu_time(t2)
c$    t2 = omp_get_wtime()

      timeinfo(1) = t2-t1
      write(*,*) timeinfo(1)
c
c     TODO: here we should probably filter out chunks which
c     are clearly worse than others for each target
c     to limit the later workload. This will increase
c     efficiency especially in the case where the search radii
c     become larger than the chunk radii
c     
      
c     convert this info (in terms of chunk per target) to
c     targets per chunk
      
      allocate(col_ptr(nch+1),row_ind(nnz),iper(nnz),iperinv(nnz),
     1     stat=ier1)
      if (ier1 .ne. 0) then
         ier = -nnz
         return
      endif
      

      call rsc_to_csc(nch,nt,nnz,row_ptr,col_ind,col_ptr,row_ind,
     1     iper)

      do i = 1,nnz
         iperinv(iper(i)) = i
      enddo
      
c     find closest chunk and point on chunk for each target

c     initialize some shared precomps
      
      k = norders(1)
      do i=2,nch
         k = max(k,norders(i))
      enddo

      ntinit = 2*k
      ncoef = k+3
      itype=0
      call legeexps(itype,ntinit,tinits,u,v,w)
      x0 = 0
      ndeg = k-1
      do i = 1,k
         pexp0(i) = 0
      enddo
      call legeexe2(x0,val0,pexp0,ndeg,pjcoefs1,pjcoefs2,ncoef)

      allocate(ts_targ2(nnz),dist22(nnz),stat=ier1)
      if (ier1 .ne. 0) then
         ier = -nnz
         return
      endif

      allocate(qwts(npts),chlens(nch))
      call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,qwts)
      do i = 1,nch
         chlens(i) = 0
         do j = ixys(i),(ixys(i+1)-1)
            chlens(i) = chlens(i)+qwts(j)
         enddo
      enddo

      
c     for now, find distance to all relevant chunks

      call cpu_time(t1)
c$    t1 = omp_get_wtime()

c$omp parallel do default(shared)
c$omp$ private(i,j,istart,iend,norderloc,xcoefs,ycoefs)
c$omp$ private(ii,jj,nt1,chlen)
      do i = 1,nch
         istart = ixys(i)
         iend = ixys(i+1)-1
         norderloc = norders(i)
         do j = 1,norderloc
            jj = istart+j-1
            xcoefs(j) = srccoefs(1,jj)
            ycoefs(j) = srccoefs(2,jj)
         enddo

         ii = col_ptr(i)
         nt1 = col_ptr(i+1)-col_ptr(i)
         ifprint = 0
         chlen = chlens(i)
         call chunknear1(norderloc,xcoefs,ycoefs,chlen,ntinit,tinits,
     1        ncoef,pjcoefs1,pjcoefs2,ndt,nt,targets,nt1,
     2        row_ind(ii),ts_targ2(ii),dist22(ii),ifprint)
      enddo
c$omp end parallel do

      
      call cpu_time(t2)
c$    t2 = omp_get_wtime()

      timeinfo(2) = t2-t1
      timeinfo(3) = nnz + 0.1d0
      
      d2max = 0
      do i=1,nnz
         d2max = max(d2max,dist22(i))
      enddo
      d2max = d2max*1.1d0 + 1.0d0
      
c     using permutation info, find closest chunk to
c     each targ and the t parameter and distance

c$omp parallel do default(shared)      
c$omp$ private(i,istart,iend,ichbest,d2best,tbest)
c$omp$ private(j,ii,d2,t)
      do i = 1,nt
         istart = row_ptr(i)
         iend = row_ptr(i+1)-1

         ichbest = -1
         d2best = d2max
         tbest = 0
         do j = istart,iend
            ii = iperinv(j)
            d2 = dist22(ii)
            t = ts_targ2(ii)
            if (d2 .lt. d2best) then
               d2best = d2
               ichbest = col_ind(j)
               tbest = t
            endif
         enddo

         ich_id(i) = ichbest
         ts_targ(i) = tbest
         dist_targ(i) = sqrt(d2best)
         if (ichbest .eq. -1) dist_targ(i) = -1
      enddo
c$omp end parallel do

      return
      end
      
      subroutine findnearchunktarg_id_ts_brute(nch,norders,ixys,iptype,
     1     npts,srccoefs,srcvals,ndt,nt,targets,ich_id,ts_targ,
     2     dist_targ)
c
c
c     Find nearest chunk and point on chunk for each target.
c     
c     This is the slower version in which all chunks are checked.
c
c     Input:
c
c     nch - integer, number of chunks
c     norders - integer(nch), discretization order of each chunk
c     ixys - integer(nch+1), index into srccoefs, srcvals arrays
c     ixys - integer(nch+1)
c       starting location in srccoefs,srcvals array
c       where info for each chunk begins
c     iptype - integer(nch)
c       type of chunk
c       iptype = 1, chunk discretized with Gauss-Legendre
c       nodes and not at a corner
c     npts - integer
c       total number of discretization points
c     srccoefs - real *8 (6,npts)
c       Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c       if(iptype.eq.1) then basis = Legendre polynomials
c     srcvals - real *8 (8,npts)
c     x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization
c       points
c     ndt - integer, first dimension of targets array
c     nt - integer, second dimension of targets array, number of targs
c     targets - real *8 (ndt,nt), array in which first two
c         entries of targets(:,i) are the coordinates of target i
c     
c     Output:
c
c     NOTE there is nothing which prevents this routine from 
c     identifying a chunk as the closest chunk to the target
c     even though the given target is in fact outside
c     of the radii specified in srcrad (in some interpolated
c     sense). I.e. sufficient closeness is not enforced.
c
c     ich_id - integer(nt), ich_id(i) is the chunk id for chunk
c         nearest to the ith target.
c         if no chunks found using the respective search radii
c         ich_id(i) = -1
c     ts_targ - real *8(nt), ts_targ(i) is the value in [-1,1]
c         of the underlying chunk parameterization of chunk ich_id(i)
c         which is closest to the target. If no close chunk was found
c         ts_targ(i) = 0
c     dist_targ - real *8(nt), dist_targ(i) is the distance
c         from the ith target to the point corresponding to 
c         ts_targ(i) on chunk ich_id(i). If no close chunk found,
c         dist_targ(i) = -1    
c      
      implicit none
      integer, intent(in) :: nch,norders(nch),ixys(nch+1),iptype(nch)
      integer, intent(in) :: npts, nt, ndt
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),
     1     targets(ndt,nt)
      real *8, intent(out) :: dist_targ(nt),ts_targ(nt)
      integer, intent(out) :: ich_id(nt)

      real *8 :: xt,yt,xs,ys,dbest,d,xcoef(100),ycoef(100),pjcoefs1(100)
      real *8 :: pjcoefs2(100),x0,val0,pexp0(100),dist21,ts_targ1,u,v
      real *8 :: tinits(200),w(200),tbest
      real *8, allocatable :: qwts(:), chlens(:)
      integer :: i,ibest,j,jj,jstart,jend,k,ntinit,ncoef,itype,ndeg
      integer :: nt1,idx(1),jjj,norderloc,ifprint
      real *8 :: chlen
c     initialize some precomps
      
      k = norders(1)
      do i=2,nch
         k = max(k,norders(i))
      enddo

      ntinit = 2*k
      ncoef = k+3
      itype=0
      call legeexps(itype,ntinit,tinits,u,v,w)
      x0 = 0
      ndeg = k-1
      do i = 1,k
         pexp0(i) = 0
      enddo
      call legeexe2(x0,val0,pexp0,ndeg,pjcoefs1,pjcoefs2,ncoef)

      allocate(qwts(npts),chlens(nch))
      call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,qwts)
      do i = 1,nch
         chlens(i) = 0
         do j = ixys(i),(ixys(i+1)-1)
            chlens(i) = chlens(i)+qwts(j)
         enddo
      enddo
      
c$omp parallel do default(shared)      
c$omp$ private(i,xt,yt,xs,ys,dbest,ibest,j,jj,jstart,jend)
c$omp$ private(nt1,idx,xcoef,ycoef,ts_targ1,dist21,norderloc)
c$omp$ private(jjj,tbest,chlen)
      do i = 1,nt
         xt = targets(1,i)
         yt = targets(2,i)

         xs = srcvals(1,1)
         ys = srcvals(2,1)         
         dbest = ((xt-xs)**2 + (yt-ys)**2)*1.0d100
         ibest = -1
         tbest = 0
         
         do j = 1,nch
            norderloc = norders(j)
            jstart = ixys(j)
            do jj = 1,norderloc
               jjj = jstart + jj -1
               xcoef(jj) = srccoefs(1,jjj)
               ycoef(jj) = srccoefs(2,jjj)
            enddo
            nt1 = 1
            idx(1) = i
            ifprint=0
            chlen = chlens(j)
            call chunknear1(norderloc,xcoef,ycoef,chlen,ntinit,tinits,
     1           ncoef,pjcoefs1,pjcoefs2,ndt,nt,targets,nt1,
     2           idx,ts_targ1,dist21,ifprint)
            if (dist21 .lt. dbest) then
               dbest = dist21
               ibest = j
               tbest = ts_targ1
            endif
         enddo

         ich_id(i) = ibest
         dist_targ(i) = sqrt(dbest)
         ts_targ(i) = tbest
      enddo
c$omp end parallel do

      return
      end


      subroutine chunknear1(norderloc,xcoefs,ycoefs,chlen,ntinit,tinits,
     1     ncoef,pjcoefs1,pjcoefs2,ndt,nt,targets,nt1,idx,ts_targ,
     2     dist2_targ,ifprint) 
c
c     find the closest point on a chunk for each target
c      
c     input
c
c     norderloc - integer, order of chunk
c     xcoefs - real *8 (norderloc), coefficients of x coordinates
c     ycoefs - real *8 (norderloc), coefficients of y coordinates      
c     chlen - real *8, length of chunk (used to scale)
c     ntinit - integer, number of elements in initial guess array
c     tinits - real *8 (ntinit), a grid in [-1,1] used
c         to find a decent starting point for algorithm. Good
c         performance has been observed using 2*norderloc
c         Legendre nodes and tacking on -1 and 1 as extras
c         (for a total ntinit = 2*norderloc + 2)
c     ncoef - integer, length of pjcoefs1,pjcoefs2 arrays,
c         must exceed norderloc
c     ndt - integer, first dimension of targets array
c     nt - integer, second dimension of targets array, number of targs
c     targets - real *8 (ndt,nt), array in which first two
c         entries of targets(:,i) are the coordinates of target i
c     nt1 - integer, number of targets of interest
c     idx - integer(nt1), indices of targets of interest in
c         targets array
c     ifprint - integer, if equal to 1 print out information
c         in the case of convergence failure 
c
c     output
c
c     ts_targ - real *8(nt1), ts_targ(i) t in [-1,1] for which
c        distance of x(t),y(t) to targets(1:2,idx(i))
c     dist2_targ - real *8(nt1), squared distance from 
c        x(ts_targ(i)),y(ts_targ(i)) to targets(1:2,idx(i))
c     
      implicit real *8 (a-h,o-z)
      real *8 :: xcoefs(norderloc),ycoefs(norderloc)

      real *8 :: pjcoefs1(ncoef),pjcoefs2(ncoef)
      real *8 :: targets(ndt,nt), tinits(ntinit)
      integer :: ncoef, norderloc, nt1, idx(nt1)
      real *8 :: ts_targ(*), dist2_targ(*)
      
c     local
      real *8 xinit(200), yinit(200), dxinit(200), dyinit(200),
     1     xcoefloc(100), ycoefloc(100)

      if (nt1 .le. 0) return
      
      ndeg = norderloc-1
      ndeg2=ndeg**2

      ifpre = 0

c     optimization parameters
      niter = 50
      iextra = 3
      thresh = 5d-8

      dlam0 = 1.0d0
      upfac = 2.0d0
      downfac = 3.0d0

      do i = 1,ntinit
         call legeexe2(tinits(i),xinit(i),xcoefs,ndeg,pjcoefs1,pjcoefs2,
     1        ifpre)
         call legeexe2(tinits(i),yinit(i),ycoefs,ndeg,pjcoefs1,pjcoefs2,
     1        ifpre)
      enddo

      do ii = 1,nt1

         i1 = idx(ii)
         
c     recentering seems to improve performance

         xt = targets(1,i1)
         yt = targets(2,i1)

c     find closest point on initial (fine-ish) grid

         x0 = xinit(1)
         y0 = yinit(1)
         d = (x0-xt)**2 + (y0-yt)**2
         dmin = d
         imin = 1

         do i = 2,ntinit
            x0 = xinit(i)
            y0 = yinit(i)
            d = (x0-xt)**2 + (y0-yt)**2
            if (d .lt. dmin) then
               dmin = d
               imin = i
            endif
         enddo

         tlow = -1
         thigh = 1
         
         t0 = tinits(imin)

c     starting from the closest, fine tune with a few iterations
c     of Levenberg algorithm

c     shift/rescale coeffs
         
         xcoefloc(1)=(xcoefs(1)-xt)/chlen
         ycoefloc(1)=(ycoefs(1)-yt)/chlen
         
         do i = 2,norderloc
            xcoefloc(i) = xcoefs(i)/chlen
            ycoefloc(i) = ycoefs(i)/chlen
         enddo

         ifend = 0

         call legefd22(t0,x0,dx0,d2x0,xcoefloc,ndeg,
     1        pjcoefs1,pjcoefs2,ifpre)
         call legefd22(t0,y0,dy0,d2y0,ycoefloc,ndeg,
     1        pjcoefs1,pjcoefs2,ifpre)

         dlam = dlam0
         
         f0 = x0**2+y0**2
         df0 = x0*dx0 + y0*dy0
         ds0 = dx0**2 + dy0**2
         d2f0 = ds0 + x0*d2x0 + y0*d2y0

c     get an upper bound for the maximum of the 
c     function xx'+yy' 

         dnrm = 0
         do i = 1,norderloc
            di2 = i*(i-1)/2.0d0
            dnrm = dnrm + (abs(xcoefloc(i))+abs(ycoefloc(i)))*(1+di2)
         enddo
         
         do i = 1,niter

c     regulate the curvature acceleration
c     we arbitrarily choose to only allow 100x the
c     levenberg step
            dreg0 = max(d2f0,1d-2*ds0)
            
            deltat1 = -df0/(dreg0*(1+dlam))
            deltat2 = -df0/(dreg0*(1+dlam/downfac))

            t1 = t0+deltat1
            t2 = t0+deltat2

            t1 = min(t1,thigh)
            t1 = max(t1,tlow)
            t2 = min(t2,thigh)
            t2 = max(t2,tlow)

            call legefd22(t1,x1,dx1,d2x1,xcoefloc,ndeg,
     1           pjcoefs1,pjcoefs2,ifpre)
            call legefd22(t1,y1,dy1,d2y1,ycoefloc,ndeg,
     1           pjcoefs1,pjcoefs2,ifpre)
            call legefd22(t2,x2,dx2,d2x2,xcoefloc,ndeg,
     1           pjcoefs1,pjcoefs2,ifpre)
            call legefd22(t2,y2,dy2,d2y2,ycoefloc,ndeg,
     1           pjcoefs1,pjcoefs2,ifpre)

            f1 = x1**2+y1**2
            f2 = x2**2+y2**2

            if (min(f1,f2) .gt. f0) then
               dlam = dlam*upfac
            else
               if (f2 .lt. f1) then
                  dlam = dlam/downfac
                  t1 = t2
                  x1 = x2
                  y1 = y2
                  dx1 = dx2
                  dy1 = dy2
                  d2x1 = d2x2
                  d2y1 = d2y2
                  f1 = f2
               endif

               t0=t1
               x0=x1
               y0=y1
               dx0=dx1
               dy0=dy1
               d2x0=d2x1
               d2y0=d2y1
               f0=f1

               f0 = x0**2+y0**2
               df0 = x0*dx0 + y0*dy0
               ds0 = dx0**2 + dy0**2
               d2f0 = ds0 + x0*d2x0 + y0*d2y0

            endif

            dkap = abs(d2f0*t0)
            denom = max(dkap,dnrm)
            der = abs(df0)
            if (der .lt. thresh*denom) ifend = ifend+1
            if (t0 .eq. tlow .and. df0 .ge. 0) ifend = iextra
            if (t0 .eq. thigh .and. df0 .le. 0) ifend = iextra 

            if (ifend .ge. iextra) exit
         enddo

         if (ifend .lt. iextra .and. ifprint .eq. 1) then
            write(*,*) "no converge"
            write(*,*) ii
            write(*,*) t0, f0, df0, d2f0
            write(*,*) der, thresh*denom
         endif

         ts_targ(ii) = t0
         dist2_targ(ii) = f0*chlen**2
      enddo

      
      return
      end
         
            

      subroutine legefd22(x,val,der,der2,pexp,n,
     1     pjcoefs1,pjcoefs2,ninit)
      implicit real *8 (a-h,o-z)
      real *8 pexp(1),pjcoefs1(1),pjcoefs2(1)
c     
C     This subroutine computes the value, derivative,
c     and second derivative of a gaussian expansion
c     with coefficients PEXP at point X in interval [-1,1].
c     
c     input parameters:
c     
C     X - evaluation point
C     PEXP - expansion coefficients
C     N  - order of expansion
c     pjcoefs1, pjcoefs2 - two arrays precomputed on a previous call
c     on a previous call to this subroutine. Please note that this
c     is only an input parameter if the parameter ninit (see below)
c     has been set to 0; otherwise, these are output parameters
c     ninit - tells the subroutine whether and to what maximum order the
c     arrays coepnm1,coepnp1,coexpnp1 should be initialized.
c     EXPLANATION: The subroutine will initialize the first ninit
c     elements of each of the arrays pjcoefs1, pjcoefs2. On the first
c     call to this subroutine, ninit should be set to the maximum
c     order n for which this subroutine might have to be called;
c     on subsequent calls, ninit should be set to 0. PLEASE NOTE
c     THAT THAT THESE ARRAYS USED BY THIS SUBROUTINE
c     ARE IDENTICAL TO THE ARRAYS WITH THE SAME NAMES USED BY THE
c     SUBROUTINE LEGEEXE2. If these arrays have been initialized
c     by one of these two subroutines, they do not need to be
c     initialized by the other one.
c     
c     IMPORTANT NOTE: n is {\bf the order of the expansion, which is
c     one less than the number of terms in the expansion!!}
c     
c     output parameters:
c     
C     VAL - computed value
C     der - computed value of the derivative
c     der2 - computed value of the second derivative derivative
C     
C     
      if(ninit .gt. 0) then
c     
        done=1
        do j=2,ninit
c     
           pjcoefs1(j)=(2*j-done)/j
           pjcoefs2(j)=-(j-done)/j
        enddo
c 
      endif

      pjm2=1
      pjm1=x
      derjm2=0
      derjm1=1
      der2jm2 = 0
      der2jm1 = 0
c     
      val=pexp(1)*pjm2+pexp(2)*pjm1
      der=pexp(2)
      der2=0
c 
      do j = 2,n
c     
         pj= pjcoefs1(j)*x*pjm1+pjcoefs2(j)*pjm2
         val=val+pexp(j+1)*pj

         derj=pjcoefs1(j)*(pjm1+x*derjm1)+pjcoefs2(j)*derjm2
         der=der+pexp(j+1)*derj

         der2j = pjcoefs1(j)*(2*derjm1+x*der2jm1)+pjcoefs2(j)*der2jm2
         der2 = der2+pexp(j+1)*der2j
         
         pjm2=pjm1
         pjm1=pj
         derjm2=derjm1
         derjm1=derj
         der2jm2=der2jm1
         der2jm1=der2j
      enddo
c 
      return
      end
c 
c 
      
