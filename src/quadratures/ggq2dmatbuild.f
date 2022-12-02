c
c     This file contains the following user callable
c     routines: 
c 
c     zchunk_matbuild_ggq - build the system matrix for the
c     supplied complex-valued integral kernel using interpolation 
c     and a generalized Gaussian quadrature (GGQ) rule
c
c     dchunk_matbuild_ggq - build the system matrix for the
c     supplied real-valued integral kernel using interpolation 
c     and a generalized Gaussian quadrature (GGQ) rule
c
c     



      subroutine zchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquad,ifrobust,
     3     sysmat,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     this subroutine builds the system matrix for the
c     supplied complex-valued integral kernel using interpolation 
c     and a generalized Gaussian quadrature (GGQ) rule
c
c  Input arguments:
c     - nch: integer
c        number of chunks
c     - norders: integer(nch)
c        order of discretization on each patch 
c     - ixys: integer(nch+1)
c        starting location of data on patch i
c     - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c     - npts: integer
c        total number of discretization points on the boundary
c     - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c     - srcvals: real *8 (8,npts)
c     x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c     - adjs - integer(2,nch)
c        adjs(1,j) is the chunk "before" chunk j 
c        adjs(2,j) is the chunk "after" chunk j
c     - fker: function handle
c         function handle for evaluating the kernel k
c         
c         expected calling sequnce:
c         fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c     - ndd: integer
c     number of real parameters for the fker routine
c     - dpars: real *8 (ndd)
c        real parameters for the fker routine
c     - ndz: integer
c         number of complex parameters for the fker routine
c     - zpars: real *8 (ndz)
c         complex parameters for the fker routine
c     - ndi: integer
c         number of integer parameters for the fker routine
c     - ipars: integer (ndi)
c         integer parameters for the fker routine
c     - opdims: integer(2)
c        the kernel is complex *16(opdims(1),opdims(2)) valued
c        for a scalar kernel opdims(1)=opdims(2)=1
c     - ising - integer
c        singularity type (ignored for now, only log implemented)
c        ising = 0, log r or weaker
c        ising = 1, 1/r weaker
c        ising = 2, 1/r^2 or weaker
c    - iquad - integer
c        quadrature type (ignored for now, only std GGQ implemented)
c        iquad = 1, use standard GGQ
c    - ifrobust - integer 
c         robustness flag (not implemented)
c         ifrobust = 1, then interactions outside of the neighbors
c                    and self which require it are treated with
c                    adaptive quadrature
c         ifrobust = 0, then this is not done
c     
c     Output arguments
c
c     sysmat - complex *16(opdims(1)*npts,opdims(2)*npts)
c              the system matrix for the supplied kernel
c     ier - integer
c        error flag
c        ier = 0, normal execution
c        ier = 1, not all chunks have same order (required for now)
c     ier = 2, requested quadrature type not available for chunk
c              orders in supplied geometry
      implicit real *8 (a-h,o-z)
      integer :: iptype,nch,norders(nch),npts,ixys(nch+1)
      real *8 :: srccoefs(6,npts), srcvals(8,npts)
      integer :: ndd, ndz, ndi, ipars(ndi), ier, adjs(2,nch)
      integer :: opdims(2), ising, iquad, ifrobust
      complex *16 :: zpars(ndz), sysmat(npts*opdims(1),npts*opdims(2))
      real *8 :: dpars(ndd)
      external :: fker
c     local
      real *8, allocatable :: xs1(:), ws1(:), xs0(:,:), ws0(:,:)
      real *8, allocatable :: ainterp1(:,:), ainterp0(:,:,:)
      real *8, allocatable :: work(:), ts(:), ws(:), dwork(:)
      integer, allocatable :: nquads0(:), iquad0lddr(:)
      complex *16, allocatable :: zwork(:), submat(:,:)
      
      ier = 0
      
      norder = norders(1)
      do i = 2,nch
         if (norders(i) .ne. norder) then
            ier = 1
            return
         endif
      enddo

c     load ggq nodes

      call ggq2dstd_get_quads_info(norder, nquad1, nquad0, ier1)
      if (ier1 .ne. 0) then
         ier = 2
         return
      endif
      allocate(xs1(nquad1), ws1(nquad1))
      allocate(xs0(nquad0,norder), ws0(nquad0,norder))
      call ggq2dstd_get_quads(norder, nquad1, xs1, ws1, nquad0,
     1     xs0, ws0)

c     this ladder structure will make it easier later
c     to use rules that have a variable number of nodes per target
      
      allocate(nquads0(norder),iquad0lddr(norder+1))
      do i = 1,norder
         nquads0(i) = nquad0
         iquad0lddr(i) = (i-1)*nquad0+1
      enddo
      iquad0lddr(norder+1) = norder*nquad0+1

c     precompute interpolation matrices

      allocate(ainterp1(nquad1,norder),ainterp0(nquad0,norder,norder))

      nquadmax = max(nquad1,nquad0)
      lwork = 2*norder**2+norder + 100
      allocate(work(lwork),ts(norder),ws(norder))
      
      call lematrin(norder,nquad1,xs1,ainterp1,ts,work)
      do i = 1,norder
         call lematrin(norder,nquad0,xs0(1,i),ainterp0(1,1,i),ts,work)
      enddo

      itype = 1
      call legeexps(itype,norder,ts,u,v,ws)

c     fill in self and adjacent entries

      idim1 = opdims(1)*norder*nch
      idim1sub = opdims(1)*norder
      idim2sub = opdims(2)*norder
      
      allocate(submat(idim1sub,idim2sub))

      nquadmax = 0
      nquadmax = max(nquadmax,nquad1)
      nquadmax = max(nquadmax,norder)
      do i = 1,norder
         nquadmax = max(nquadmax,nquads0(i))
      enddo
      
      allocate(dwork((8+1)*nquadmax),
     1     zwork(opdims(1)*opdims(2)*nquadmax*norder))

      ndt = 8
      ntarg = norder

      do jj = 1,nch
         iileft = adjs(1,jj)
         iiright = adjs(1,jj)
         jstart = (jj-1)*norder + 1
         jstartmat = (jj-1)*norder*opdims(2) + 1
         do ii = 1,nch
            istart = (ii-1)*norder + 1
            istartmat = (ii-1)*norder*opdims(1) + 1
            
            if (ii .eq. jj) then
               call zchunk_indrulematbuild_ggq(norder,srcvals(1,jstart),
     1              ndt,ntarg,srcvals(1,jstart),nquads0,
     2              iquad0lddr,ainterp0,ws0,fker,ndd,dpars,ndz,zpars,
     3              ndi,ipars,opdims,submat,zwork,dwork)
               call zchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else if (ii .eq. iileft) then
               call zchunk_onerulematbuild_ggq(norder,srcvals(1,jstart),
     1              nquad1,ainterp1,ws1,ndt,ntarg,srcvals(1,istart),
     2              fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              zwork,dwork)
               call zchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else if (ii .eq. iiright) then
               call zchunk_onerulematbuild_ggq(norder,srcvals(1,jstart),
     1              nquad1,ainterp1,ws1,ndt,ntarg,srcvals(1,istart),
     2              fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              zwork,dwork)
               call zchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else
               call zchunk_nativematbuild(norder,srcvals(1,jstart),
     1              ws,ndt,ntarg,srcvals(1,istart),fker,
     2              ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              zwork,dwork)
               call zchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            endif
         enddo
      enddo
      
      return
      end

      subroutine zchunk_buildmat_writein(nrow,a,nrow1,ncol1,a1)
      implicit none
      integer nrow,nrow1,ncol1,i,j
      complex *16 :: a(nrow,*), a1(nrow1,ncol1)

      do j = 1,ncol1
         do i = 1,nrow1
            a(i,j) = a1(i,j)
         enddo
      enddo
      
      return
      end

      subroutine chunk_interpsrcvals(norder,srcvals,
     1     nquad,ainterp,srcvalsinterp,dsdtinterp)
      implicit real *8 (a-h,o-z)
      real *8 :: srcvals(8,norder), srcvalsinterp(8,nquad),
     1     dsdtinterp(nquad)
      real *8 :: ainterp(nquad,norder)

      do j = 1,nquad
         srcvalsinterp(1,j) = 0
         srcvalsinterp(2,j) = 0
         srcvalsinterp(3,j) = 0
         srcvalsinterp(4,j) = 0
         srcvalsinterp(5,j) = 0
         srcvalsinterp(6,j) = 0
         srcvalsinterp(7,j) = 0
         srcvalsinterp(8,j) = 0
      enddo
      
      do i = 1,norder
         sc1 = srcvals(1,i)
         sc2 = srcvals(2,i)
         sc3 = srcvals(3,i)
         sc4 = srcvals(4,i)
         sc5 = srcvals(5,i)
         sc6 = srcvals(6,i)
         sc7 = srcvals(7,i)
         sc8 = srcvals(8,i)         
         do j = 1,nquad
            srcvalsinterp(1,j) = srcvalsinterp(1,j) + ainterp(j,i)*sc1
            srcvalsinterp(2,j) = srcvalsinterp(2,j) + ainterp(j,i)*sc2
            srcvalsinterp(3,j) = srcvalsinterp(3,j) + ainterp(j,i)*sc3
            srcvalsinterp(4,j) = srcvalsinterp(4,j) + ainterp(j,i)*sc4
            srcvalsinterp(5,j) = srcvalsinterp(5,j) + ainterp(j,i)*sc5
            srcvalsinterp(6,j) = srcvalsinterp(6,j) + ainterp(j,i)*sc6
            srcvalsinterp(7,j) = srcvalsinterp(7,j) + ainterp(j,i)*sc7
            srcvalsinterp(8,j) = srcvalsinterp(8,j) + ainterp(j,i)*sc8
         enddo
      enddo

      do j = 1,nquad
         dx = srcvalsinterp(3,j)
         dy = srcvalsinterp(4,j)
         ds = sqrt(dx**2+dy**2)
         srcvalsinterp(7,j) = dy/ds
         srcvalsinterp(8,j) = -dx/ds
         dsdtinterp(j) = ds
      enddo
      
      return
      end


      subroutine zchunk_applyrule_onetarg_ggq(ndt,targinfo,norder,
     1     nquad,srcinterp,dsdtinterp,ws,ainterp,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,nrowsub,submat,zval)
      implicit real *8 (a-h,o-z)
      real *8 targinfo(ndt), srcinterp(8,nquad), ws(nquad),
     1     ainterp(nquad,norder), dsdtinterp(nquad)
      external fker
      integer ndd,ndz,ndi
      real *8 dpars(*)
      integer ipars(*), opdims(2)
      complex *16 :: zpars(*), zval(opdims(1),*)
      complex *16 :: submat(nrowsub,*)

      do ii = 1,nquad
         ioff = (ii-1)*opdims(2)         
         call fker(srcinterp(1,ii),ndt,targinfo,ndd,dpars,
     1        ndz,zpars,ndi,ipars,zval(1,ioff+1))
         do j = 1,opdims(2)
            do i = 1,opdims(1)
               zval(i,ioff+j) = zval(i,ioff+j)*ws(ii)*dsdtinterp(ii)
            enddo
         enddo
      enddo
      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         do j = 1,opdims(2)
            do i = 1,opdims(1)
               submat(i,joff+j)=0
            enddo
         enddo
      enddo

      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         do ii = 1,nquad
            ioff = (ii-1)*opdims(2)
            da = ainterp(ii,jj)
            do j = 1,opdims(2)
               do i = 1,opdims(1)
                  submat(i,joff+j) = submat(i,joff+j)
     1                 + zval(i,ioff+j)*da
               enddo
            enddo
         enddo
      enddo

      return
      end
      
      subroutine zchunk_indrulematbuild_ggq(norder,srcvals,
     1     ndt,ntarg,targs,nquads0,iquad0lddr,ainterp0,ws0,
     2     fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3     zwork,dwork)
c
c     build submatrix with an individual rule for each
c     target (typically this is for the self interaction)
c
      implicit real *8 (a-h,o-z)
      external fker
      real *8 :: srcvals(8,*), ws0(*), ainterp0(*), dwork(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder, nquads0(ntarg), iquad0lddr(ntarg+1)
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      complex *16 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zwork(*), zpars(*)

      nrowsub = opdims(1)*ntarg
      maxquad = 0
      do i = 1,ntarg
         maxquad = max(maxquad,nquads0(i))
      enddo
      lsrcinterp = 8*maxquad
      isrcinterp = 1
      idsdtinterp = isrcinterp+lsrcinterp

      do ii = 1,ntarg
         iquad = iquad0lddr(ii)
         nquad = nquads0(ii)
         iinterp = (iquad-1)*norder + 1
         isub = (ii-1)*opdims(1)+1
         call chunk_interpsrcvals(norder,srcvals,
     1        nquad,ainterp0(iinterp),dwork(isrcinterp),
     2        dwork(idsdtinterp))
         call zchunk_applyrule_onetarg_ggq(ndt,targs(1,ii),
     1        norder,nquad,dwork(isrcinterp),dwork(idsdtinterp),
     2        ws0(iquad),ainterp0(iinterp),fker,ndd,dpars,ndz,zpars,
     2        ndi,ipars,opdims,nrowsub,submat(isub,1),zwork)

      enddo

      return
      end

      subroutine zchunk_onerulematbuild_ggq(norder,srcvals,
     1     nquad,ainterp,ws,ndt,ntarg,targs,fker,
     2     ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3     zwork,dwork)
c
c     build submatrix using the same rule for all
c     targets (typically this is for neighbors or
c     over-sampled interactions)
c
      implicit real *8 (a-h,o-z)
      external fker
      real *8 :: srcvals(8,*), ws(*), dwork(*), ainterp(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      complex *16 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zwork(*),zpars(*)

      nrowsub = opdims(1)*ntarg
      lsrcinterp = 8*nquad
      isrcinterp = 1
      idsdtinterp = isrcinterp+lsrcinterp

      call chunk_interpsrcvals(norder,srcvals,nquad,ainterp,
     1     dwork(isrcinterp),dwork(idsdtinterp))
      
      do ii = 1,ntarg
         isub = (ii-1)*opdims(1)+1
         call zchunk_applyrule_onetarg_ggq(ndt,targs(1,ii),
     1        norder,nquad,dwork(isrcinterp),dwork(idsdtinterp),
     2        ws,ainterp,fker,ndd,dpars,ndz,zpars,
     2        ndi,ipars,opdims,nrowsub,submat(isub,1),zwork)
      enddo

      
      return
      end

      subroutine zchunk_nativematbuild(norder,srcvals,
     1     ws,ndt,ntarg,targs,fker,ndd,dpars,ndz,zpars,
     2     ndi,ipars,opdims,submat,zval,dsdt)
c
c     build submatrix using the native rule for all
c     targets (similar to one rule case but no interpolation
c     needed)
c
      implicit real *8 (a-h,o-z)
      real *8 :: srcvals(8,*), ws(*), dsdt(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      complex *16 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zval(opdims(1),*), zpars(*)

      do i = 1,norder
         dsdt(i) = sqrt(srcvals(3,i)**2 + srcvals(4,i)**2)
      enddo

      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         dw = ws(jj)
         ds = dsdt(jj)
         do ii = 1,ntarg
            ioff = (ii-1)*opdims(1)
            call fker(srcvals(1,jj),ndt,targs(1,ii),ndd,dpars,
     1           ndz,zpars,ndi,ipars,zval)
            do j = 1,opdims(2)
               do i = 1,opdims(1)
                  submat(ioff+i,joff+j) = zval(i,j)*dw*ds
               enddo
            enddo
         enddo
      enddo

      return
      end

      
      subroutine dchunk_matbuild_ggq(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,ising,iquad,ifrobust,
     3     sysmat,ier)
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c     this subroutine builds the system matrix for the
c     supplied real-valued integral kernel using interpolation 
c     and a generalized Gaussian quadrature (GGQ) rule
c
c  Input arguments:
c     - nch: integer
c        number of chunks
c     - norders: integer(nch)
c        order of discretization on each patch 
c     - ixys: integer(nch+1)
c        starting location of data on patch i
c     - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c     - npts: integer
c        total number of discretization points on the boundary
c     - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c     - srcvals: real *8 (8,npts)
c     x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c     - adjs - integer(2,nch)
c        adjs(1,j) is the chunk "before" chunk j 
c        adjs(2,j) is the chunk "after" chunk j
c     - fker: function handle
c         function handle for evaluating the kernel k
c         
c         expected calling sequnce:
c         fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c     - ndd: integer
c     number of real parameters for the fker routine
c     - dpars: real *8 (ndd)
c        real parameters for the fker routine
c     - ndz: integer
c         number of complex parameters for the fker routine
c     - zpars: real *8 (ndz)
c         complex parameters for the fker routine
c     - ndi: integer
c         number of integer parameters for the fker routine
c     - ipars: integer (ndi)
c         integer parameters for the fker routine
c     - opdims: integer(2)
c        the kernel is real *8(opdims(1),opdims(2)) valued
c        for a scalar kernel opdims(1)=opdims(2)=1
c     - ising - integer
c        singularity type (ignored for now, only log implemented)
c        ising = 0, log r or weaker
c        ising = 1, 1/r weaker
c        ising = 2, 1/r^2 or weaker
c    - iquad - integer
c        quadrature type (ignored for now, only std GGQ implemented)
c        iquad = 1, use standard GGQ
c    - ifrobust - integer 
c         robustness flag (not implemented)
c         ifrobust = 1, then interactions outside of the neighbors
c                    and self which require it are treated with
c                    adaptive quadrature
c         ifrobust = 0, then this is not done
c     
c     Output arguments
c
c     sysmat - complex *16(opdims(1)*npts,opdims(2)*npts)
c              the system matrix for the supplied kernel
c     ier - integer
c        error flag
c        ier = 0, normal execution
c        ier = 1, not all chunks have same order (required for now)
c     ier = 2, requested quadrature type not available for chunk
c              orders in supplied geometry
c      
      implicit real *8 (a-h,o-z)
      integer :: iptype,nch,norders(nch),npts,ixys(nch+1)
      real *8 :: srccoefs(6,npts), srcvals(8,npts)
      integer :: ndd, ndz, ndi, ipars(ndi), ier, adjs(2,nch)
      integer :: opdims(2), ising, iquad, ifrobust
      real *8 :: sysmat(npts*opdims(1),npts*opdims(2))
      complex *16 :: zpars(ndz)
      real *8 :: dpars(ndd)
      external :: fker
c     local
      real *8, allocatable :: xs1(:), ws1(:), xs0(:,:), ws0(:,:)
      real *8, allocatable :: ainterp1(:,:), ainterp0(:,:,:)
      real *8, allocatable :: work(:), ts(:), ws(:), dwork(:)
      integer, allocatable :: nquads0(:), iquad0lddr(:)
      real *8, allocatable :: submat(:,:)
      
      
      ier = 0
      
      norder = norders(1)
      do i = 2,nch
         if (norders(i) .ne. norder) then
            ier = 1
            return
         endif
      enddo

c     load ggq nodes

      call ggq2dstd_get_quads_info(norder, nquad1, nquad0,ier1)
      if (ier1 .ne. 0) then
         ier = 2
         return
      endif
      
      allocate(xs1(nquad1), ws1(nquad1))
      allocate(xs0(nquad0,norder), ws0(nquad0,norder))
      call ggq2dstd_get_quads(norder, nquad1, xs1, ws1, nquad0,
     1     xs0, ws0)

c     this ladder structure will make it easier later
c     to use rules that have a variable number of nodes per target
      
      allocate(nquads0(norder),iquad0lddr(norder+1))
      do i = 1,norder
         nquads0(i) = nquad0
         iquad0lddr(i) = (i-1)*nquad0+1
      enddo
      iquad0lddr(norder+1) = norder*nquad0+1

c     precompute interpolation matrices

      allocate(ainterp1(nquad1,norder),ainterp0(nquad0,norder,norder))

      nquadmax = max(nquad1,nquad0)
      lwork = 2*norder**2+norder + 100
      allocate(work(lwork),ts(norder),ws(norder))
      
      call lematrin(norder,nquad1,xs1,ainterp1,ts,work)
      do i = 1,norder
         call lematrin(norder,nquad0,xs0(1,i),ainterp0(1,1,i),ts,work)
      enddo

      itype = 1
      call legeexps(itype,norder,ts,u,v,ws)

c     fill in self and adjacent entries

      idim1 = opdims(1)*norder*nch
      idim1sub = opdims(1)*norder
      idim2sub = opdims(2)*norder
      
      allocate(submat(idim1sub,idim2sub))

      nquadmax = 0
      nquadmax = max(nquadmax,nquad1)
      nquadmax = max(nquadmax,norder)
      do i = 1,norder
         nquadmax = max(nquadmax,nquads0(i))
      enddo

      lwork= (8+1)*nquadmax + opdims(1)*opdims(2)*nquadmax*norder
      allocate(dwork(lwork))

      idsdt = 1
      ldsdt = nquadmax
      ikernvals = idsdt+ldsdt
      
      ndt = 8
      ntarg = norder

      do jj = 1,nch
         iileft = adjs(1,jj)
         iiright = adjs(1,jj)
         jstart = (jj-1)*norder + 1
         jstartmat = (jj-1)*norder*opdims(2) + 1
         do ii = 1,nch
            istart = (ii-1)*norder + 1
            istartmat = (ii-1)*norder*opdims(1) + 1
            
            if (ii .eq. jj) then
               call dchunk_indrulematbuild_ggq(norder,srcvals(1,jstart),
     1              ndt,ntarg,srcvals(1,jstart),nquads0,
     2              iquad0lddr,ainterp0,ws0,fker,ndd,dpars,ndz,zpars,
     3              ndi,ipars,opdims,submat,dwork)
               call dchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else if (ii .eq. iileft) then
               call dchunk_onerulematbuild_ggq(norder,srcvals(1,jstart),
     1              nquad1,ainterp1,ws1,ndt,ntarg,srcvals(1,istart),
     2              fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              dwork)
               call dchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else if (ii .eq. iiright) then
               call dchunk_onerulematbuild_ggq(norder,srcvals(1,jstart),
     1              nquad1,ainterp1,ws1,ndt,ntarg,srcvals(1,istart),
     2              fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              dwork)
               call dchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            else
               call dchunk_nativematbuild(norder,srcvals(1,jstart),
     1              ws,ndt,ntarg,srcvals(1,istart),fker,
     2              ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3              dwork(ikernvals),dwork)
               call dchunk_buildmat_writein(idim1,
     1              sysmat(istartmat,jstartmat),idim1sub,idim2sub,
     2              submat)
            endif
         enddo
      enddo
      
      return
      end

      subroutine dchunk_buildmat_writein(nrow,a,nrow1,ncol1,a1)
      implicit none
      integer nrow,nrow1,ncol1,i,j
      real *8 :: a(nrow,*), a1(nrow1,ncol1)

      do j = 1,ncol1
         do i = 1,nrow1
            a(i,j) = a1(i,j)
         enddo
      enddo
      
      return
      end


      subroutine dchunk_applyrule_onetarg_ggq(ndt,targinfo,norder,
     1     nquad,srcinterp,dsdtinterp,ws,ainterp,fker,ndd,dpars,
     2     ndz,zpars,ndi,ipars,opdims,nrowsub,submat,dval)
      implicit real *8 (a-h,o-z)
      real *8 targinfo(ndt), srcinterp(8,nquad), ws(nquad),
     1     ainterp(nquad,norder), dsdtinterp(nquad)
      external fker
      integer ndd,ndz,ndi
      real *8 dpars(*)
      integer ipars(*), opdims(2)
      real *8 :: dval(opdims(1),*)
      complex *16 :: zpars(*)
      real *8 :: submat(nrowsub,*)

      do ii = 1,nquad
         ioff = (ii-1)*opdims(2)         
         call fker(srcinterp(1,ii),ndt,targinfo,ndd,dpars,
     1        ndz,zpars,ndi,ipars,dval(1,ioff+1))
         do j = 1,opdims(2)
            do i = 1,opdims(1)
               dval(i,ioff+j) = dval(i,ioff+j)*ws(ii)*dsdtinterp(ii)
            enddo
         enddo
      enddo
      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         do j = 1,opdims(2)
            do i = 1,opdims(1)
               submat(i,joff+j)=0
            enddo
         enddo
      enddo

      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         do ii = 1,nquad
            ioff = (ii-1)*opdims(2)
            da = ainterp(ii,jj)
            do j = 1,opdims(2)
               do i = 1,opdims(1)
                  submat(i,joff+j) = submat(i,joff+j)
     1                 + dval(i,ioff+j)*da
               enddo
            enddo
         enddo
      enddo

      return
      end
      
      subroutine dchunk_indrulematbuild_ggq(norder,srcvals,
     1     ndt,ntarg,targs,nquads0,iquad0lddr,ainterp0,ws0,
     2     fker,ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3     dwork)
c
c     build submatrix with an individual rule for each
c     target (typically this is for the self interaction)
c
      implicit real *8 (a-h,o-z)
      external fker
      real *8 :: srcvals(8,*), ws0(*), ainterp0(*), dwork(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder, nquads0(ntarg), iquad0lddr(ntarg+1)
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      real *8 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zpars(*)

      nrowsub = opdims(1)*ntarg
      maxquad = 0
      do i = 1,ntarg
         maxquad = max(maxquad,nquads0(i))
      enddo
      lsrcinterp = 8*maxquad
      isrcinterp = 1
      idsdtinterp = isrcinterp+lsrcinterp
      ldsdtinterp = maxquad
      ikernvals = idsdtinterp + ldsdtinterp
      
      do ii = 1,ntarg
         iquad = iquad0lddr(ii)
         nquad = nquads0(ii)
         iinterp = (iquad-1)*norder + 1
         isub = (ii-1)*opdims(1)+1
         call chunk_interpsrcvals(norder,srcvals,
     1        nquad,ainterp0(iinterp),dwork(isrcinterp),
     2        dwork(idsdtinterp))
         call dchunk_applyrule_onetarg_ggq(ndt,targs(1,ii),
     1        norder,nquad,dwork(isrcinterp),dwork(idsdtinterp),
     2        ws0(iquad),ainterp0(iinterp),fker,ndd,dpars,ndz,zpars,
     2        ndi,ipars,opdims,nrowsub,submat(isub,1),
     4        dwork(ikernvals))

      enddo

      return
      end

      subroutine dchunk_onerulematbuild_ggq(norder,srcvals,
     1     nquad,ainterp,ws,ndt,ntarg,targs,fker,
     2     ndd,dpars,ndz,zpars,ndi,ipars,opdims,submat,
     3     dwork)
c
c     build submatrix using the same rule for all
c     targets (typically this is for neighbors or
c     over-sampled interactions)
c
      implicit real *8 (a-h,o-z)
      external fker
      real *8 :: srcvals(8,*), ws(*), dwork(*), ainterp(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      real *8 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zpars(*)

      nrowsub = opdims(1)*ntarg
      lsrcinterp = 8*nquad
      isrcinterp = 1
      idsdtinterp = isrcinterp+lsrcinterp
      ldsdtinterp = nquad
      ikernvals = idsdtinterp+ldsdtinterp

      call chunk_interpsrcvals(norder,srcvals,nquad,ainterp,
     1     dwork(isrcinterp),dwork(idsdtinterp))
      
      do ii = 1,ntarg
         isub = (ii-1)*opdims(1)+1
         call dchunk_applyrule_onetarg_ggq(ndt,targs(1,ii),
     1        norder,nquad,dwork(isrcinterp),dwork(idsdtinterp),
     2        ws,ainterp,fker,ndd,dpars,ndz,zpars,
     2        ndi,ipars,opdims,nrowsub,submat(isub,1),
     3        dwork(ikernvals))
      enddo

      
      return
      end

      subroutine dchunk_nativematbuild(norder,srcvals,
     1     ws,ndt,ntarg,targs,fker,ndd,dpars,ndz,zpars,
     2     ndi,ipars,opdims,submat,dval,dsdt)
c
c     build submatrix using the native rule for all
c     targets (similar to one rule case but no interpolation
c     needed)
c
      implicit real *8 (a-h,o-z)
      real *8 :: srcvals(8,*), ws(*), dsdt(*)
      real *8 :: targs(ndt,ntarg), dpars(*)
      integer :: norder
      integer :: ndd,ndz,ndi,ipars(*),opdims(2)
      real *8 :: submat(opdims(1)*ntarg,opdims(2)*norder)
      complex *16 :: zpars(*)
      real *8 :: dval(opdims(1),*)

      do i = 1,norder
         dsdt(i) = sqrt(srcvals(3,i)**2 + srcvals(4,i)**2)
      enddo

      do jj = 1,norder
         joff = (jj-1)*opdims(2)
         dw = ws(jj)
         ds = dsdt(jj)
         do ii = 1,ntarg
            ioff = (ii-1)*opdims(1)
            call fker(srcvals(1,jj),ndt,targs(1,ii),ndd,dpars,
     1           ndz,zpars,ndi,ipars,dval)
            do j = 1,opdims(2)
               do i = 1,opdims(1)
                  submat(ioff+i,joff+j) = dval(i,j)*dw*ds
               enddo
            enddo
         enddo
      enddo

      return
      end
      



      subroutine zgetoncurvequad_ggq2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,eps,ipv,
     2     fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     3     col_ind,iquad,nquad,wnear,ier)
c
c------------------------
c     This subroutine generates the near field quadrature
c     for on curve interactions
c     for a given kernel which is assumed to be
c     a compact value integral operator (principal value and 
c     hypersingular options will be added in later)
c     where the near field is specified by the user 
c     in row sparse compressed format.
c     
c     The quadrature is computed by generalized Gaussian
c     quadrature for self and near interactions and adaptive
c     integration for other interactions
c
c  Input arguments:
c
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each chunk
c    - ixys: integer(nch+1)
c        ixys(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to chunk i
c    - iptype: integer(nch)
c        type of chunk
c        *  iptype = 1, triangular chunk discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (6,npts)
c        Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = Legendre polynomials
c    - srcvals: double precision (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at discretization points
c    - eps: double precision
c        precision requested
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators (Currently not supported)
c        * ipv = 2, for hypersingular operators (currently not supported)
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nnz: integer
c        number of source chunk-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source chunkes
c        for target i start
c    - col_ind: integer (nnz)
c        list of source chunkes relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c
c  Output parameters:
c
c    - wnear: double complex(nquad)
c        near field quadrature corrections
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndi,ndd,ndz,ipv
      integer, intent(out) :: ier, adjs(2,nch)
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: nch,norders(nch),npts
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer nchmax
      integer, allocatable :: row_ind(:),col_ptr(:),iper(:),ichunk(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),ts(:),wts(:)
      complex *16, allocatable :: zints(:),ztmp(:)
      integer itarg,kmax,opdims(2)
      real *8, allocatable :: xs1(:), ws1(:), xs0(:,:), ws0(:,:)
      real *8, allocatable :: ainterp1(:,:), ainterp0(:,:,:)
      real *8, allocatable :: work(:), dwork(:),srcnbor(:,:)
      integer, allocatable :: nquads0(:), iquad0lddr(:), inborlist(:)
      complex *16, allocatable :: zwork(:), submat(:,:),submatnbor(:)
      
      external fker

      ier = 0
      
      norder = norders(1)
      do i = 2,nch
         if (norders(i) .ne. norder) then
            ier = 1
            return
         endif
      enddo

      itype = 2
      allocate(ts(norder),umat(norder,norder),vmat(norder,norder),
     1     wts(norder))
      call legeexps(itype,norder,ts,umat,vmat,wts)
      allocate(zints(norder+5),ztmp(norder+5))

c     load ggq nodes

      call ggq2dstd_get_quads_info(norder, nquad1, nquad0, ier1)
      if (ier1 .ne. 0) then
         ier = 2
         return
      endif
      allocate(xs1(nquad1), ws1(nquad1))
      allocate(xs0(nquad0,norder), ws0(nquad0,norder))
      call ggq2dstd_get_quads(norder, nquad1, xs1, ws1, nquad0,
     1     xs0, ws0)

c     this ladder structure will make it easier later
c     to use rules that have a variable number of nodes per target
      
      allocate(nquads0(norder),iquad0lddr(norder+1))
      do i = 1,norder
         nquads0(i) = nquad0
         iquad0lddr(i) = (i-1)*nquad0+1
      enddo
      iquad0lddr(norder+1) = norder*nquad0+1

c     precompute interpolation matrices

      allocate(ainterp1(nquad1,norder),ainterp0(nquad0,norder,norder))

      nquadmax = max(nquad1,nquad0)
      lwork = 2*norder**2+norder + 100
      allocate(work(lwork))
      
      call lematrin(norder,nquad1,xs1,ainterp1,ts,work)
      do i = 1,norder
         call lematrin(norder,nquad0,xs0(1,i),ainterp0(1,1,i),ts,work)
      enddo

c     get transpose of interaction info

      allocate(row_ind(nnz),col_ptr(nch+1),iper(nnz))
      call rsc_to_csc(nch,npts,nnz,row_ptr,col_ind,col_ptr,
     1     row_ind,iper)
      
C$OMP PARALLEL DO
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO      

      allocate(ichunk(npts))
      
c$OMP PARALLEL DO
      do i = 1,nch
         do j = ixys(i),(ixys(i+1)-1)
            ichunk(j)=i
         enddo
      enddo
c$OMP END PARALLEL DO      

c     storage for ggq work

      opdims(1)=1
      opdims(2)=1
      
      idim1 = opdims(1)*norder*nch
      idim1sub = opdims(1)*norder
      idim1nbor = opdims(1)*norder*2
      idim2sub = opdims(2)*norder
      
      allocate(submat(idim1sub,idim2sub))

      nquadmax = 0
      nquadmax = max(nquadmax,nquad1)
      nquadmax = max(nquadmax,norder)
      do i = 1,norder
         nquadmax = max(nquadmax,nquads0(i))
      enddo
      
      allocate(dwork((8+1)*nquadmax),
     1     zwork(opdims(1)*opdims(2)*nquadmax*norder),
     2     inborlist(2*norder),srcnbor(8,2*norder),
     3     submatnbor(idim1nbor*idim2sub))


c     adaptive integration params
      
      nqorder = 10
      nchmax = 10000
      ndt = 8

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itarg,j,ich,istart,zints)
C$OMP$PRIVATE(ztmp,ntarg0,zwork,dwork,submat,ichtarg,ichl,ichr,irel)
C$OMP$PRIVATE(k,inborlist,inbor,nnbor,srcnbor,submatnbor)
      do ich=1,nch

         ichl=adjs(1,ich)
         ichr=adjs(2,ich)
         istart = ixys(ich)
         ntarg0=norder
         call zchunk_indrulematbuild_ggq(norder,srcvals(1,istart),
     1        ndt,ntarg0,srcvals(1,istart),nquads0,
     2        iquad0lddr,ainterp0,ws0,fker,ndd,dpars,ndz,zpars,
     3        ndi,ipars,opdims,submat,zwork,dwork)
         
         nnbor=0
         do j=col_ptr(ich),col_ptr(ich+1)-1
            itarg = row_ind(j)
            ichtarg = ichunk(itarg)
            if (ichtarg .eq. ichl .or. ichtarg .eq. ichr) then
               nnbor=nnbor+1
               inborlist(nnbor)=itarg
               do k=1,8
                  srcnbor(k,nnbor)=srcvals(k,itarg)
               enddo
            endif
         enddo

         if (nnbor .gt. 0) then
            call zchunk_onerulematbuild_ggq(norder,srcvals(1,istart),
     1           nquad1,ainterp1,ws1,ndt,nnbor,srcnbor,fker,
     2           ndd,dpars,ndz,zpars,ndi,ipars,opdims,submatnbor,
     3           zwork,dwork)
         endif
         
         inbor=0
         do j=col_ptr(ich),col_ptr(ich+1)-1
            itarg = row_ind(j)
            ichtarg = ichunk(itarg)
            if (ichtarg .eq. ich) then
c     self
               irel = itarg-istart+1
               do k=1,norder
                  wnear(iquad(iper(j))+k-1)=submat(irel,k)
               enddo
            else if (ichtarg .eq. ichl .or. ichtarg .eq. ichr) then
c     neighbor rule
               inbor=inbor+1
               call zchunk_buildmat_readrow(nnbor,norder,
     1              submatnbor,inbor,wnear(iquad(iper(j))))
            else
               ntarg0=1
               call zchunkints_lege(eps,norders(ich),srccoefs(1,istart),
     1              ndt,ntarg0,srcvals(1,itarg),norders(ich),fker,
     2              ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nchmax,
     3              zints)
               call zrmatmatt_f77(1,norders(ich),zints,
     1              norders(ich),umat,wnear(iquad(iper(j))))
               
            endif
         enddo
      enddo
C$OMP END PARALLEL DO      


c
c
      return
      end
c
c
c
c
      
      subroutine zchunk_buildmat_readrow(nrow,ncol,a,irow,r)
      implicit none
      integer nrow,ncol,irow,j
      complex *16 :: a(nrow,ncol), r(ncol)

      do j = 1,ncol
         r(j) = a(irow,j)
      enddo
      
      return
      end


      subroutine dgetoncurvequad_ggq2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,eps,ipv,
     2     fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,
     3     col_ind,iquad,nquad,wnear,ier)
c
c------------------------
c     This subroutine generates the near field quadrature
c     for on curve interactions
c     for a given kernel which is assumed to be
c     a compact value integral operator (principal value and 
c     hypersingular options will be added in later)
c     where the near field is specified by the user 
c     in row sparse compressed format.
c     
c     The quadrature is computed by generalized Gaussian
c     quadrature for self and near interactions and adaptive
c     integration for other interactions
c
c  Input arguments:
c
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each chunk
c    - ixys: integer(nch+1)
c        ixys(i) denotes the starting location in srccoefs,
c        and srcvals array corresponding to chunk i
c    - iptype: integer(nch)
c        type of chunk
c        *  iptype = 1, triangular chunk discretized using RV nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (6,npts)
c        Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = Legendre polynomials
c    - srcvals: double precision (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at discretization points
c    - eps: double precision
c        precision requested
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators (Currently not supported)
c        * ipv = 2, for hypersingular operators (currently not supported)
c    - fker: procedure pointer
c        function handle for the kernel. Calling sequence 
c        * fker(srcinfo,ndtarg,targinfo,ndd,dpars,ndz,zpars,ndi,ipars,val)
c    - ndd: integer
c        number of double precision parameters
c    - dpars: double precision(ndd)
c        double precision parameters
c    - ndz: integer
c        number of double complex parameters
c    - zpars: double complex(ndz)
c        double complex parameters
c    - ndi: integer
c        number of integer parameters
c    - ipars: integer(ndi)
c        integer parameters
c    - nnz: integer
c        number of source chunk-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer
c        to col_ind array where list of relevant source chunkes
c        for target i start
c    - col_ind: integer (nnz)
c        list of source chunkes relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c
c  Output parameters:
c
c    - wnear: double (nquad)
c        near field quadrature corrections
c----------------------------------------------------               
c
      implicit real *8 (a-h,o-z)
      integer, intent(in) :: ndi,ndd,ndz,ipv
      integer, intent(out) :: ier, adjs(2,nch)
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: nch,norders(nch),npts
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      real *8, intent(out) :: wnear(nquad)

      integer nchmax
      integer, allocatable :: row_ind(:),col_ptr(:),iper(:),ichunk(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),ts(:),wts(:)
      real *8, allocatable :: dints(:),dtmp(:)
      integer itarg,kmax,opdims(2)
      real *8, allocatable :: xs1(:), ws1(:), xs0(:,:), ws0(:,:)
      real *8, allocatable :: ainterp1(:,:), ainterp0(:,:,:)
      real *8, allocatable :: work(:), dwork(:), srcnbor(:,:)
      integer, allocatable :: nquads0(:), iquad0lddr(:), inborlist(:)
      real *8, allocatable :: submat(:,:),submatnbor(:)
      
      external fker

      ier = 0
      
      norder = norders(1)
      do i = 2,nch
         if (norders(i) .ne. norder) then
            ier = 1
            return
         endif
      enddo

      itype = 2
      allocate(ts(norder),umat(norder,norder),vmat(norder,norder),
     1     wts(norder))
      call legeexps(itype,norder,ts,umat,vmat,wts)
      allocate(dints(norder+5),dtmp(norder+5))

c     load ggq nodes

      call ggq2dstd_get_quads_info(norder, nquad1, nquad0, ier1)
      if (ier1 .ne. 0) then
         ier = 2
         return
      endif
      allocate(xs1(nquad1), ws1(nquad1))
      allocate(xs0(nquad0,norder), ws0(nquad0,norder))
      call ggq2dstd_get_quads(norder, nquad1, xs1, ws1, nquad0,
     1     xs0, ws0)

c     this ladder structure will make it easier later
c     to use rules that have a variable number of nodes per target
      
      allocate(nquads0(norder),iquad0lddr(norder+1))
      do i = 1,norder
         nquads0(i) = nquad0
         iquad0lddr(i) = (i-1)*nquad0+1
      enddo
      iquad0lddr(norder+1) = norder*nquad0+1

c     precompute interpolation matrices

      allocate(ainterp1(nquad1,norder),ainterp0(nquad0,norder,norder))

      nquadmax = max(nquad1,nquad0)
      lwork = 2*norder**2+norder + 100
      allocate(work(lwork))
      
      call lematrin(norder,nquad1,xs1,ainterp1,ts,work)
      do i = 1,norder
         call lematrin(norder,nquad0,xs0(1,i),ainterp0(1,1,i),ts,work)
      enddo

c     get transpose of interaction info
      
      allocate(row_ind(nnz),col_ptr(nch+1),iper(nnz))
      call rsc_to_csc(nch,npts,nnz,row_ptr,col_ind,col_ptr,
     1     row_ind,iper)

C$OMP PARALLEL DO
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO      

      allocate(ichunk(npts))
      
c$OMP PARALLEL DO
      do i = 1,nch
         do j = ixys(i),(ixys(i+1)-1)
            ichunk(j)=i
         enddo
      enddo
c$OMP END PARALLEL DO      

c     storage for ggq work

      opdims(1)=1
      opdims(2)=1
      
      idim1 = opdims(1)*norder*nch
      idim1sub = opdims(1)*norder
      idim1nbor = opdims(1)*norder*2
      idim2sub = opdims(2)*norder
      
      allocate(submat(idim1sub,idim2sub))

      nquadmax = 0
      nquadmax = max(nquadmax,nquad1)
      nquadmax = max(nquadmax,norder)
      do i = 1,norder
         nquadmax = max(nquadmax,nquads0(i))
      enddo

      lwork= (8+1)*nquadmax + opdims(1)*opdims(2)*nquadmax*norder      
      allocate(dwork(lwork),
     2     inborlist(2*norder),srcnbor(8,2*norder),
     3     submatnbor(idim1nbor*idim2sub))


c     adaptive integration params
      
      nqorder = 10
      nchmax = 10000

c     always 8 because these are source points
      ndt = 8

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itarg,j,ich,istart,dints)
C$OMP$PRIVATE(dtmp,ntarg0,dwork,submat,ichtarg,ichl,ichr,irel)
C$OMP$PRIVATE(k,inborlist,inbor,nnbor,srcnbor,submatnbor)
      do ich=1,nch

         ichl=adjs(1,ich)
         ichr=adjs(2,ich)
         istart = ixys(ich)
         ntarg0=norder
         call dchunk_indrulematbuild_ggq(norder,srcvals(1,istart),
     1        ndt,ntarg0,srcvals(1,istart),nquads0,
     2        iquad0lddr,ainterp0,ws0,fker,ndd,dpars,ndz,zpars,
     3        ndi,ipars,opdims,submat,dwork)
         
         nnbor=0
         do j=col_ptr(ich),col_ptr(ich+1)-1
            itarg = row_ind(j)
            ichtarg = ichunk(itarg)
            if (ichtarg .eq. ichl .or. ichtarg .eq. ichr) then
               nnbor=nnbor+1
               inborlist(nnbor)=itarg
               do k=1,8
                  srcnbor(k,nnbor)=srcvals(k,itarg)
               enddo
            endif
         enddo

         if (nnbor .gt. 0) then
            call dchunk_onerulematbuild_ggq(norder,srcvals(1,istart),
     1           nquad1,ainterp1,ws1,ndt,nnbor,srcnbor,fker,
     2           ndd,dpars,ndz,zpars,ndi,ipars,opdims,submatnbor,
     3           dwork)
         endif

         inbor=0
         do j=col_ptr(ich),col_ptr(ich+1)-1
            itarg = row_ind(j)
            ichtarg = ichunk(itarg)
            if (ichtarg .eq. ich) then
c     self
               irel = itarg-istart+1
               do k=1,norder
                  wnear(iquad(iper(j))+k-1)=submat(irel,k)
               enddo
            else if (ichtarg .eq. ichl .or. ichtarg .eq. ichr) then
c     neighbor rule
               inbor=inbor+1
               call dchunk_buildmat_readrow(nnbor,norder,
     1              submatnbor,inbor,wnear(iquad(iper(j))))
            else
c     any other target by adaptive
               ntarg0=1
               call dchunkints_lege(eps,norders(ich),srccoefs(1,istart),
     1              ndt,ntarg0,srcvals(1,itarg),norders(ich),fker,
     2              ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nchmax,
     3              dints)
               call dmatmatt(1,norders(ich),dints,
     1              norders(ich),umat,wnear(iquad(iper(j))))
              
            endif
         enddo
      enddo
C$OMP END PARALLEL DO      


c
c
      return
      end
c
c
      
      subroutine dchunk_buildmat_readrow(nrow,ncol,a,irow,r)
      implicit none
      integer nrow,ncol,irow,j
      real *8 :: a(nrow,ncol), r(ncol)

      do j = 1,ncol
         r(j) = a(irow,j)
      enddo
      
      return
      end

      
