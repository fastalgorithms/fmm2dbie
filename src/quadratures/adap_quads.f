c
c
c     This file has the following user callable routines
c
c        ?getnearquad_ggq_guru - guru interface for 
c         computing near field quadrature for a compact/pv kernel
c        
c
c         z - complex
c
c
c
c
c
c

      subroutine zgetnearquad_adap_guru2d(nch,norders,ixys,iptype,npts,
     1   srccoefs,srcvals,ndtarg,ntarg,targvals,ich_id,ts_targ,eps,ipv,
     2   fker,ndd,dpars,ndz,zpars,ndi,ipars,nnz,row_ptr,col_ind,iquad,
     3   nquad,wnear)
c
c------------------------
c  This subroutine generates the near field quadrature
c  for a given kernel which is assumed to be
c  a compact value integral operator (principal value and 
c  hypersingular options will be added in later)
c  where the near field is specified by the user 
c  in row sparse compressed format.
c
c
c  The quadrature is computed by adaptive integration
c  for targets off the chunk and adaptive integration
c  by splitting into 2 chunks.
c  
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
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targvals: double precision (ndtarg,ntarg)
c        target info. First two components must be x,y
c        coordinates
c    - ich_id: integer(ntarg)
c        chunk id of target, ich_id = -1, if target off-surface
c    - ts_targ: double precision (ntarg)
c        local coordinates on chunk if target on surface
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
      integer, intent(in) :: ipars(ndi)
      integer, intent(in) :: ndtarg
      integer, intent(in) :: nch,norders(nch),npts
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      integer, intent(in) :: ntarg,ich_id(ntarg)
      real *8, intent(in) :: ts_targ(ntarg)
      real *8, intent(in) :: targvals(ndtarg,ntarg)
      real *8, intent(in) :: dpars(ndd)
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)

      integer nchmax
      integer, allocatable :: row_ind(:),col_ptr(:),iper(:)
      real *8, allocatable :: umat(:,:),vmat(:,:),ts(:),wts(:)
      complex *16, allocatable :: zints(:)
      integer itarg,kmax

      external fker

      kmax = maxval(norders)
      allocate(zints(kmax+5))

      allocate(row_ind(nnz),col_ptr(nch+1),iper(nnz))
      call rsc_to_csc(nch,ntarg,nnz,row_ptr,col_ind,col_ptr,row_ind,
     1  iper)

C$OMP PARALLEL DO
      do i=1,quad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO      

      ntarg0 = 1

      nqorder = 10
      nchmax = 10000
      itype = 2

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(itarg,j,ich,istart,zints)
C$OMP$PRIVATE(ts,umat,vmat,wts,kk)
      do ich=1,nch
        kk = norders(ich)
        allocate(ts(kk),umat(kk,kk),vmat(kk,kk),wts(kk))
        call legeexps(itype,kk,ts,umat,vmat,wts)
        
        istart = ixys(ich)
        do j=col_ptr(ich),col_ptr(ich+1)-1
          itarg = row_ind(j)

          if(iptype(ich).eq.1) then
            if(ich.ne.ich_id(itarg)) then
              call zchunkints_lege(eps,norders(ich),srccoefs(1,istart),
     1          ndtarg,ntarg0,targvals(1,itarg),norders(ich),fker,
     2          ndd,dpars,ndz,zpars,ndi,ipars,nqorder,nchmax,zints)
            else
              call zget_adap_lege_self(ipv,eps,norders(ich), 
     1          ts_targ(itarg),srccoefs(1,istart),ndtarg,
     1          targvals(1,itarg),norders(ich),fker,ndd,
     2          dpars,ndz,zpars,ndi,ipars,nqorder,zints)
            endif
          endif
          call zrmatmatt(1,norders(ich),zints,norders(ich),umat,
     1       wnear(iquad(j)))
        enddo
        deallocate(ts,umat,vmat,wts)
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
c
      subroutine zget_adap_lege_self(ipv,eps,norder,tt0,srccoefs,ndtarg,
     1   xytarg,nporder,fker,ndd,dpars,ndz,zpars,ndi,
     2   ipars,nchmax,zintvals)
c
c   This subrourine computes the integrals
c      
c     \int_{-1}^{1} K(x_{i},\rho(t)) P_{n}(t) dsdt(t) dt \, ,
c  
c   where $P_{n}(t)$ are the legendre polynomials,
c   $\rho(t)$ is the parametrization of a chunk
c   dsdt is the speed of parametrization of $\rho(t)$,
c   and $x_{i}$ are the target locations.
c
c   The map $\rho(t)$ is stored as coefficients of a legendre expansion
c   $xy(t)$ along with expansions of the first and second derivative
c   information with respect to t.
c
c   The integrals are computed using adaptive integration 
c   by splitting the chunk into 2 intervals  until an absolute 
c   tolerance of $\eps$ is reached.
c
c   Input arguments:
c    - ipv: integer
c        Flag for choosing type of self-quadrature
c        * ipv = 0, for compact/weakly singular operators
c        * ipv = 1, for singular operators (Currently not supported)
c        * ipv = 2, for hypersingular operators (currently not supported)
c     - eps: real *8
c         requested tolerance
c     - tt0: real *8 
c         location of target on the panel
c     - norder: integer
c         order of discretization of input chunk
c     - srccoefs: real *8 (6,norder)
c         Legendre expansion coefficients of x,y,dxdt,dydt,d2xdt2,d2ydt2
c     - ndtarg: integer
c         leading order dimension of target location array
c     - xytarg: real *8 (ndtarg)
c         target info. The first two components xytarg(1:2,i) must be 
c         the xy coordinates of the targets
c     - nporder: integer
c         order of legendre polynomials to be integrated
c     - fker: function handle
c         function handle for evaluating the kernel k
c         
c         expected calling sequnce:
c         fker(x,ndtarg,y,ndd,dpars,ndz,zpars,ndi,ipars,f)
c     - ndd: integer
c         number of real parameters for the fker routine
c     - dpars: real *8 (ndd)
c         real parameters for the fker routine
c     - ndz: integer
c         number of complex parameters for the fker routine
c     - zpars: real *8 (ndz)
c         complex parameters for the fker routine
c     - ndi: integer
c         number of integer parameters for the fker routine
c     - ipars: integer (ndi)
c         integer parameters for the fker routine
c     - nqorder: integer 
c         order of quadrature to be used for adaptive integration
c     - nchmax: integer
c         maximum number of chunks in adaptive integration, 
c         will be increased by factors of 2 after max is reached
c         for current chunk (Feature under construction)
c    
c   Output arguments:
c     - zintvals: complex *16 (nporder)
c         computed integrals
c   
c-----------------------------------
      implicit none
      integer, intent(in) :: ipv
      real *8, intent(in) :: eps,tt0
      integer, intent(in) :: norder,ndtarg,nporder,ndd,ndz,ndi
      integer, intent(in) :: ipars(ndi),nchmax
      real *8, intent(in) :: srccoefs(6,norder),dpars(ndd)
      real *8, intent(in) :: xytarg(ndtarg)
      complex *16, intent(in) :: zpars(ndz)

      complex *16, intent(out) :: zintvals(nporder)

      real *8, allocatable :: stack(:,:)
      complex *16, allocatable :: vals(:,:)
      complex *16, allocatable :: value2(:),value3(:)
      complex *16, allocatable :: ztmp(:)
      real *8, allocatable :: ts0(:),w0(:)
      real *8 umat0,vmat0
      integer itype,nqorder

      real *8 a,b
      integer maxrec,numit,maxdepth,nnmax,i,j,l,ier
      integer iord

      external fker

      itype = 1

      if(norder.lt.5) then
        iord = 1
      else if(norder.lt.10) then
        iord = 2
      else if(norder.lt.15) then
        iord = 3
      else
        iord = 4
      endif

      iord = 1


      call get_nqorder_selfquad(ipv,iord,nqorder)
      call prinf('nqorder=*',nqorder,1)

      allocate(ts0(nqorder),w0(nqorder))
      call load_selfquad(ipv,iord,nqorder,ts0,w0)
      call prin2('ts0=*',ts0,nqorder)
      call prin2('w0=*',w0,nqorder)

      maxdepth = 200
      allocate(stack(2,maxdepth),vals(nporder,maxdepth))
      allocate(value2(nporder),value3(nporder))
      allocate(ztmp(nporder))

      a = -1.0d0
      b = tt0

      nnmax = nchmax
        
      ier = 0
      maxrec = 0
      numit = 0
      call zadinrecm(ier,stack,a,b,norder,srccoefs,ndtarg,xytarg,
     1  nporder,fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ts0,w0,
     2     vals,nnmax,eps,zintvals,maxdepth,maxrec,numit,value2,
     3     value3)
      print *, "numit=",numit

      a = tt0
      b = 1.0d0
      ier = 0
      maxrec = 0
      numit = 0
      call zadinrecm(ier,stack,a,b,norder,srccoefs,ndtarg,xytarg,
     1  nporder,fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ts0,w0,
     2     vals,nnmax,eps,ztmp,maxdepth,maxrec,numit,value2,
     3     value3)
      print *, "numit=",numit
      
      do i=1,nporder
        zintvals(i) = zintvals(i) + ztmp(i)
      enddo

      return
      end
c
c
c
