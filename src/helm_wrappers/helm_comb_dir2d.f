
c
c     This file contains the following user callable
c     routines: 
c 
c       getnearquad_helm_comb_dir_2d - generates the near
c        field quadrature for the Dirichlet data
c        corresponding to the combined field
c        representation 
c
c       lpcomp_helm_comb_dir_2d 
c          simpler version of helmholtz layer potential evaluator
c          only geometry, targets, representation parameters (alpha,beta,k)
c          and density sampled at discretization required on input,
c          output is the layer potential evaluated at the target points
c          (note that the identity term is not included for targets on
c           surface)
c
c       helm_comb_dir_solver_2d - solves the interior/exterior Dirichlet
c         problem for Helmholtz equation using the combined field
c         representation
c
c
c    Advanced user interfaces: 
c       lpcomp_helm_comb_dir_addsub_2d 
c         compute layer potential for the Dirichlet
c         data corresponding to the combined field representation
c         using add and subtract
c 
c
c
c



      subroutine getnearquad_helm_comb_dir_2d(nch,norders,
     1   ixys,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ich_id,ts_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3   iquad,nquad,wnear)
c
c------------------------------
c  This subroutine generates the near field quadrature 
c  for the representation 
c
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the quadrature.
c
c  Currently, the only scheme available for computing the quadrature
c  is adaptive integration
c
c  Input arguments:
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: real *8 (ndtarg,ntarg)
c        target information, the first two components must be
c        xy coordinates
c    - ich_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - ts_targ: real *8 (ntarg)
c        local t coordinate of chunk if on-surface
c    - eps: real *8
c        precision requested
c    - zpars: complex *16 (3)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
c        zpars(3) = beta
c    - iquadtype - integer
c        quadrature type
c        iquadtype = 1, use adaptive quadrature for everything 
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer to col_ind array where list of 
c        relevant source patches for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c
c  Output
c
c    - wnear: complex *16(nquad)
c        the desired near field quadrature
c               
c

      implicit none 
      integer, intent(in) :: nch,norders(nch),npts,nquad
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: iquadtype
      real *8, intent(in) :: targs(ndtarg,ntarg)
      integer, intent(in) :: ich_id(ntarg)
      real *8, intent(in) :: ts_targ(ntarg)
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(ntarg+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)


      integer ipars
      integer ndd,ndz,ndi
      real *8 dpars

      complex *16 alpha,beta
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external h2d_slp, h2d_dlp, h2d_comb

c
c
c        initialize the appropriate kernel function
c

      alpha = zpars(2)
      beta = zpars(3)

      ndz = 3
      ndi = 0
      ndd = 0
      if(iquadtype.eq.1) then
        fker => h2d_comb
        if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
          fker=>h2d_slp
        else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
          fker=>h2d_dlp
        endif
        ipv = 0
        

        call zgetnearquad_adap_guru2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ich_id,ts_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     1     ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,wnear)
      endif


      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*alpha
        enddo
C$OMP END PARALLEL DO        
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*beta
        enddo
C$OMP END PARALLEL DO        
      endif

      return
      end
c
c
c
c
c
      subroutine lpcomp_helm_comb_dir_2d(nch,norders,ixys,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   ich_id,ts_targ,eps,zpars,sigma,pot)
c
cf2py intent(in) nch,norders,ixys,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ich_id,ts_targ,eps,zpars
cf2py intent(in) sigma
cf2py intent(out) pot
c
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation 
c
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the layer potential.
c
c
c  Input arguments:
c
c    - nch: integer
c        number of patches
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        ixys(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer(nch)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (6,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (8,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: double precision (ndtarg,ntarg)
c        target information
c    - ich_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - ts_targ: double precision (2,ntarg)
c        local uv coordinates on patch if on surface, otherwise
c        set to 0 by default
c    - eps: double precision
c        precision requested
c    - zpars: double complex (3)
c        kernel parameters (Referring to formula above)
c        zpars(1) = k 
c        zpars(2) = $\alpha$
c        zpars(3) = $\beta$
c     - sigma: double complex(npts)
c         density for layer potential
c
c  Output arguments
c    - pot: double complex(ntarg)
c        layer potential evaluated at the target points
c
c-----------------------------------
c
      implicit none
      integer, intent(in) :: nch,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(nch),ixys(nch+1)
      integer, intent(in) :: iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(3)
      complex *16, intent(in) :: sigma(npts)
      integer, intent(in) :: ich_id(ntarg)
      real *8, intent(in) :: ts_targ(ntarg)

      complex *16, intent(out) :: pot(ntarg)


      integer nptso,nnz,nquad


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 over4pi
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
      data over4pi/0.07957747154594767d0/


      iptype_avg = floor(sum(iptype)/(nch+0.0d0))
      norder_avg = floor(sum(norders)/(nch+0.0d0))

      call get_rfac2d(norder_avg,iptype_avg,rfac) 
      allocate(cms(2,nch),rads(nch),rad_near(nch))

      call get_centroid_rads2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,nch
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c
      call findnear2dmem(cms,nch,rad_near,ndtarg,targs,ntarg,nnz)

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear2d(cms,nch,rad_near,ndtarg,targs,ntarg,row_ptr, 
     1        col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc2d(nch,ixys,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

      ikerorder = -1
      if(abs(dpars(2)).gt.1.0d-16) ikerorder = 0

c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(nch),ixyso(nch+1))

c      call get_far_order2d(eps,nch,norders,ixys,iptype,cms,
c     1    rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),
c     2    nnz,row_ptr,col_ind,rfac,novers,ixyso)

      do i=1,nch
        novers(i) = norders(i)
        ixyso(i) = ixys(i)
      enddo
      ixyso(nch+1) = ixys(nch+1)
      npts_over = ixyso(nch+1)-1

      allocate(srcover(8,npts_over),wover(npts_over))

      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyso,npts_over,srcover)

      call get_qwts2d(nch,novers,ixyso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call getnearquad_helm_comb_dir_2d(nch,norders,
     1      ixys,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1      ich_id,ts_targ,eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1      iquad,nquad,wnear)


c
c
c   compute layer potential
c
      call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1  iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2  eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3  sigma,novers,npts_over,ixyso,srcover,wover,pot)



      return
      end
c
c
c
c
c

      subroutine lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1   iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     2   eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,sigma,novers,
     3   nptso,ixyso,srcover,whtsover,pot)
c
cf2py intent(in) nch,norders,ixys,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,eps,zpars
cf2py intent(in) sigma,nnz,row_ptr,col_ind,iquad,nquad,wnear,novers
cf2py intent(in) nptso,ixyso,srcover,whtsover
cf2py intent(out) pot
c
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation 
c
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  This is a guru interface for the layer potential evaluator.
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the layer potential.
c
c
c  Input arguments:
c
c    - nch: integer
c        number of patches
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        ixys(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer(nch)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (6,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (8,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: double precision (ndtarg,ntarg)
c        target information
c    - eps: double precision
c        precision requested
c    - zpars: double complex (3)
c        kernel parameters (Referring to formula above)
c        zpars(1) = k 
c        zpars(2) = $\alpha$
c        zpars(3) = $\beta$
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer to col_ind array where list of 
c        relevant source patches for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c    - wnear: complex *16(nquad)
c        the desired near field quadrature
c    - sigma: double complex(npts)
c        density for layer potential
c    - novers: integer(nch)
c        order of discretization for oversampled sources and
c        density
c    - nptso: integer
c        number of oversampled points
c    - ixyso: integer(nch+1)
c        ixyso(i) denotes the starting location in srcover,
c        whtsover array, and other functions sampled at oversampled
c        nodes where information for patch i begins
c    - srcover: double precision (8,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        oversampled nodes
c    - whtsover: double precision (npts)
c        quadrature weights for integrating smooth functions sampled
c        at the oversampled nodes
c
c  Output arguments
c    - pot: double complex(ntarg)
c        layer potential evaluated at the target points
c
c-----------------------------------
c               
c
      implicit none
      integer, intent(in) :: nch,npts
      integer, intent(in) :: ndtarg,ntarg
      integer, intent(in) :: norders(nch),ixys(nch+1)
      integer, intent(in) :: ixyso(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      real *8, intent(in) :: targs(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz,row_ptr(ntarg+1),col_ind(nnz),nquad
      integer, intent(in) :: iquad(nnz+1)
      complex *16, intent(in) :: wnear(nquad),sigma(npts)
      integer, intent(in) :: novers(nch)
      integer, intent(in) :: nptso
      real *8, intent(in) :: srcover(8,nptso),whtsover(nptso)
      complex *16, intent(out) :: pot(ntarg)

      integer norder,npols,nover,npolso
      complex *16, allocatable :: potsort(:)

      real *8, allocatable :: sources(:,:),targvals(:,:)
      complex *16, allocatable :: charges(:),dipstr(:),sigmaover(:)
      real *8, allocatable :: dipvec(:,:)
      integer ns,nt
      complex *16 alpha,beta
      integer ifcharge,ifdipole
      integer ifpgh,ifpghtarg
      complex *16 tmp(10),val

      real *8 xmin,xmax,ymin,ymax,zmin,zmax,sizey,sizez,boxsize


      integer i,j,jpatch,jquadstart,jstart


      integer ifaddsub

      integer ntj
      
      complex *16 zdotu,pottmp
      real *8 radexp,epsfmm

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime

      real *8, allocatable :: radsrc(:)
      real *8, allocatable :: srctmp2(:,:)
      complex *16, allocatable :: ctmp2(:),dtmp2(:)
      real *8, allocatable :: dipvec2(:,:)
      real *8 thresh,ra
      real *8 rr,rmin
      real *8 over4pi
      integer nss,ii,l,npover
      integer nmax,ier,iper

      integer nd,ntarg0

      real *8 ttot,done,pi

      parameter (nd=1,ntarg0=1)

      ns = nptso
      done = 1
      pi = atan(done)*4

c
c    estimate max number of sources in neear field of 
c    any target
c
      nmax = 0
      call get_near_corr_max2d(ntarg,row_ptr,nnz,col_ind,nch,
     1  ixyso,nmax)
      allocate(srctmp2(2,nmax),ctmp2(nmax),dtmp2(nmax))
      allocate(dipvec2(2,nmax))
           
      ifpgh = 0
      ifpghtarg = 1
      allocate(sources(2,ns),targvals(2,ntarg))
      allocate(charges(ns),dipstr(ns),dipvec(2,ns))
      allocate(sigmaover(ns))

c 
c       oversample density
c
      
      call oversample_fun_curv2d(2,nch,norders,ixys,iptype, 
     1    npts,sigma,novers,ixyso,ns,sigmaover)

      ra = 0
c
c       set relevatn parameters for the fmm
c
      alpha = zpars(2)
      beta = zpars(3)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)      
      do i=1,ns
        sources(1,i) = srcover(1,i)
        sources(2,i) = srcover(2,i)

        charges(i) = sigmaover(i)*whtsover(i)*alpha
        dipstr(i) = sigmaover(i)*whtsover(i)*beta
        dipvec(1,i) = srcover(7,i)
        dipvec(2,i) = srcover(8,i)
      enddo
C$OMP END PARALLEL DO      

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        targvals(1,i) = targs(1,i)
        targvals(2,i) = targs(2,i)
      enddo
C$OMP END PARALLEL DO      

      ifcharge = 1
      ifdipole = 1

      if(alpha.eq.0) ifcharge = 0
      if(beta.eq.0) ifdipole = 0

c
c
c       call the fmm
c
      call cpu_time(t1)
C$      t1 = omp_get_wtime()      
      call hfmm2d(nd,eps,zpars(1),ns,sources,ifcharge,charges,
     1  ifdipole,dipstr,dipvec,iper,ifpgh,tmp,tmp,tmp,ntarg,
     1  targvals,ifpghtarg,pot,tmp,tmp,ier)
      call cpu_time(t2)
C$      t2 = omp_get_wtime()

            
      timeinfo(1) = t2-t1


c
c        compute threshold for ignoring local computation
c
      call get_fmm2d_thresh(2,ns,sources,2,ntarg,targvals,thresh)


c
c       add in precomputed quadrature

      call cpu_time(t1)
C$      t1 = omp_get_wtime()

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,jquadstart)
C$OMP$PRIVATE(jstart,pottmp,npols,l)
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          npols = ixys(jpatch+1)-ixys(jpatch)
          jquadstart = iquad(j)
          jstart = ixys(jpatch) 
          do l=1,npols
            pot(i) = pot(i) + wnear(jquadstart+l-1)*sigma(jstart+l-1)
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO

c

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,jpatch,srctmp2)
C$OMP$PRIVATE(ctmp2,dtmp2,dipvec2,nss,l,jstart,ii,val,npover)
      do i=1,ntarg
        nss = 0
        do j=row_ptr(i),row_ptr(i+1)-1
          jpatch = col_ind(j)
          do l=ixyso(jpatch),ixyso(jpatch+1)-1
            nss = nss+1
            srctmp2(1,nss) = srcover(1,l)
            srctmp2(2,nss) = srcover(2,l)

            if(ifcharge.eq.1) ctmp2(nss) = charges(l)
            if(ifdipole.eq.1) then
              dtmp2(nss) = dipstr(l)
              dipvec2(1,nss) = dipvec(1,l)
              dipvec2(2,nss) = dipvec(2,l)
            endif
          enddo
        enddo

        val = 0
        if(ifcharge.eq.1.and.ifdipole.eq.0) then
          call h2d_directcp(nd,zpars(1),srctmp2,nss,ctmp2,
     1        targvals(1,i),1,val,thresh)
        endif

        if(ifcharge.eq.0.and.ifdipole.eq.1) then
          call h2d_directdp(nd,zpars(1),srctmp2,nss,dtmp2,
     1          dipvec2,targvals(1,i),1,val,thresh)
        endif

        if(ifcharge.eq.1.and.ifdipole.eq.1) then
          call h2d_directcdp(nd,zpars(1),srctmp2,nss,ctmp2,dtmp2,
     1          dipvec2,targvals(1,i),1,val,thresh)
        endif
        pot(i) = pot(i) - val
      enddo
      
      call cpu_time(t2)
C$      t2 = omp_get_wtime()     

      timeinfo(2) = t2-t1


c      call prin2('lpcomp timeinfo=*',timeinfo,2)
c      call prinf('nfmm *',ns,1)
c      call prinf('nquad *',nquad,1)      
      ttot = timeinfo(1) + timeinfo(2)
cc      call prin2('time in lpcomp=*',ttot,1)

      
      return
      end

c
c
c
c
c
c
c
c        
      subroutine helm_comb_dir_solver_2d(nch,norders,ixys,
     1    iptype,npts,srccoefs,srcvals,adjs,eps,zpars,numit,ifinout,
     2    rhs,eps_gmres,niter,errs,rres,soln)
c
c  Solve the Helmholtz boundary value problem using the combined 
c  field integral equation
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c
c  Input arguments:
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - eps: real *8
c        precision requested
c    - zpars: complex *16 (3)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
c        zpars(3) = beta
c    - numit: integer
c        max number of gmres iterations
c    - ifinout: integer
c        flag for interior or exterior problems (normals assumed to 
c        be pointing in exterior of region)
c        * ifinout = 0, interior problem
c        * ifinout = 1, exterior problem
c    - rhs: complex *16(npts)
c        right hand side
c    - eps_gmres: real *8
c        gmres tolerance requested
c 
c  Output arguments:
c    - niter: integer
c        number of gmres iterations required for relative residual
c    - errs: real *8 (1:niter)
c        relative residual as a function of iteration number
c    - rres: real *8
c        relative residual for computed solution
c    - soln: complex *16(npts)
c        density which solves the dirichlet problem
c-----------------------------------
c
      implicit none
      integer nch,norder,npols,npts
      integer ifinout
      integer norders(nch),ixys(nch+1)
      integer iptype(nch), adjs(2,nch)
      real *8 srccoefs(6,npts),srcvals(8,npts),eps,eps_gmres
      complex *16 zpars(3)
      complex *16 rhs(npts)
      complex *16 soln(npts)

      real *8, allocatable :: targs(:,:)
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_targ(:)
      integer ndtarg,ntarg

      real *8 errs(numit+1)
      real *8 rres,eps2
      integer niter


      integer nover,npolso,nptso
      integer nnz,nquad
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8 :: rho
      integer :: ier, ising, npolyfac
      real *8, allocatable :: srcover(:,:),wover(:),rects(:,:,:)
      integer, allocatable :: ixyso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0,t0
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over

c
c
c       gmres variables
c
      complex *16 zid,ztmp
      real *8 rb,wnrm2
      integer numit,it,iind,it1,k,l
      real *8 rmyerr
      complex *16 temp
      complex *16, allocatable :: vmat(:,:),hmat(:,:)
      complex *16, allocatable :: cs(:),sn(:)
      complex *16, allocatable :: svec(:),yvec(:),wtmp(:)


      allocate(vmat(npts,numit+1),hmat(numit,numit))
      allocate(cs(numit),sn(numit))
      allocate(wtmp(npts),svec(numit+1),yvec(numit+1))


      done = 1
      pi = atan(done)*4

c
c        setup targets as on surface discretization points
c 
      ndtarg = 2
      ntarg = npts
      allocate(targs(ndtarg,npts))

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,ntarg
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
      enddo
C$OMP END PARALLEL DO   


      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0

c
c     define near field and oversampling
c

      rho = 2d0
      npolyfac=2

      allocate(novers(nch),ixyso(nch+1))
      
      call ellipse_nearfield2d_getnovers(eps,rho,npolyfac,
     1     nch,norders,ising,novers,ier)

      ixyso(1)=1
      do i = 1,nch
         ixyso(i+1)=ixyso(i)+novers(i)
      enddo
      npts_over=ixyso(nch+1)-1

      allocate(rects(2,4,nch))
      call ellipse_nearfield2d_definerects(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,rho,rects)

      call findinrectangle_mem(nch,rects,npts,ndtarg,targs,nnz,
     1     ier)

      allocate(row_ptr(npts+1),col_ind(nnz))

      call findinrectangle(nch,rects,npts,ndtarg,targs,row_ptr,
     1     nnz,col_ind,ier)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,
     1     iquad)

c
c     oversample geometry and get oversampled weights
c     for FMM
c      
      
      allocate(srcover(8,npts_over),wover(npts_over))

      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyso,npts_over,srcover)

      call get_qwts2d(nch,novers,ixyso,iptype,npts_over,
     1        srcover,wover)

c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call cpu_time(t0)
      
      call getoncurvequad_helm_comb_dir_2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,
     2     eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3     iquad,nquad,wnear)
      

      call cpu_time(t1)

      call prin2('time generating near quad *',t1-t0,1)
      
      print *, "done generating near quadrature, now starting gmres"


c
c
c     start gmres code here
c
c     NOTE: matrix equation should be of the form (z*I + K)x = y
c       the identity scaling (z) is defined via zid below,
c       and K represents the action of the principal value 
c       part of the matvec
c
      zid = -(-1)**(ifinout)*zpars(3)/2


      niter=0

c
c      compute norm of right hand side and initialize v
c 
      rb = 0

      do i=1,numit
        cs(i) = 0
        sn(i) = 0
      enddo


c
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rb)
      do i=1,npts
        rb = rb + abs(rhs(i))**2
      enddo
C$OMP END PARALLEL DO      
      rb = sqrt(rb)

C$OMP PARALLEL DO DEFAULT(SHARED)
      do i=1,npts
        vmat(i,1) = rhs(i)/rb
      enddo
C$OMP END PARALLEL DO      

      svec(1) = rb

      do it=1,numit
        it1 = it + 1

c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


        call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1       iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2       eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3       vmat(1,it),novers,npts_over,ixyso,srcover,wover,wtmp)

        do k=1,it
          ztmp = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:ztmp)          
          do j=1,npts
            ztmp = ztmp + wtmp(j)*conjg(vmat(j,k))
          enddo
C$OMP END PARALLEL DO          
          hmat(k,it) = ztmp

C$OMP PARALLEL DO DEFAULT(SHARED) 
          do j=1,npts
            wtmp(j) = wtmp(j)-hmat(k,it)*vmat(j,k)
          enddo
C$OMP END PARALLEL DO          
        enddo
          
        hmat(it,it) = hmat(it,it)+zid
        wnrm2 = 0
C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:wnrm2)        
        do j=1,npts
          wnrm2 = wnrm2 + abs(wtmp(j))**2
        enddo
C$OMP END PARALLEL DO        
        wnrm2 = sqrt(wnrm2)

C$OMP PARALLEL DO DEFAULT(SHARED) 
        do j=1,npts
          vmat(j,it1) = wtmp(j)/wnrm2
        enddo
C$OMP END PARALLEL DO        

        do k=1,it-1
          temp = cs(k)*hmat(k,it)+conjg(sn(k))*hmat(k+1,it)
          hmat(k+1,it) = -sn(k)*hmat(k,it)+cs(k)*hmat(k+1,it)
          hmat(k,it) = temp
        enddo

        ztmp = wnrm2

        call zrotmat_gmres2d(hmat(it,it),ztmp,cs(it),sn(it))
          
        hmat(it,it) = cs(it)*hmat(it,it)+conjg(sn(it))*wnrm2
        svec(it1) = -sn(it)*svec(it)
        svec(it) = cs(it)*svec(it)
        rmyerr = abs(svec(it1))/rb
        errs(it) = rmyerr
        print *, "iter=",it,errs(it)

        if(rmyerr.le.eps_gmres.or.it.eq.numit) then

c
c            solve the linear system corresponding to
c            upper triangular part of hmat to obtain yvec
c
c            y = triu(H(1:it,1:it))\s(1:it);
c
          do j=1,it
            iind = it-j+1
            yvec(iind) = svec(iind)
            do l=iind+1,it
              yvec(iind) = yvec(iind) - hmat(iind,l)*yvec(l)
            enddo
            yvec(iind) = yvec(iind)/hmat(iind,iind)
          enddo



c
c          estimate x
c
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j)
          do j=1,npts
            soln(j) = 0
            do i=1,it
              soln(j) = soln(j) + yvec(i)*vmat(j,i)
            enddo
          enddo
C$OMP END PARALLEL DO          


          rres = 0
C$OMP PARALLEL DO DEFAULT(SHARED)          
          do i=1,npts
            wtmp(i) = 0
          enddo
C$OMP END PARALLEL DO          
c
c        NOTE:
c        replace this routine by appropriate layer potential
c        evaluation routine  
c


          call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1      iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2      eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3      soln,novers,npts_over,ixyso,srcover,wover,wtmp)

C$OMP PARALLEL DO DEFAULT(SHARED) REDUCTION(+:rres)            
          do i=1,npts
            rres = rres + abs(zid*soln(i) + wtmp(i)-rhs(i))**2
          enddo
C$OMP END PARALLEL DO          
          rres = sqrt(rres)/rb
          niter = it
          return
        endif
      enddo
c
      return
      end
c
c
c
c
c
c        



      subroutine getnearquad_helm_comb_dir_2d_matlab(nch,norders,
     1   ixys,iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,cms,
     2   rad_near,nnz,nquad,zk,wnear,iind,jind)
c
c------------------------------
c  This subroutine generates the near field quadrature 
c  for the representation 
c
c
c  .. math ::
c  
c      u = (i \mathcal{S}_{k} + \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the quadrature.
c
c  Currently, the only scheme available for computing the quadrature
c  is adaptive integration
c
c  Input arguments:
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: real *8 (ndtarg,ntarg)
c        target information, the first two components must be
c        xy coordinates
c    - nquad: integer
c        number of entries in wnear
c
c  Output
c
c    - wnear: complex *16(nquad)
c        the desired near field quadrature
c    - iind: integer(nquad)
c        row indices of wnear
c    - jind: integer(nquad)
c        column indices of wnear
c               
c

      implicit none 
      integer, intent(in) :: nch,norders(nch),npts,nquad
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts)
      integer, intent(in) :: ndtarg,ntarg
      real *8, intent(in) :: targs(ndtarg,ntarg)
      real *8, intent(in) :: cms(2,nch),rad_near(nch)
      integer, intent(in) :: nnz
      integer, allocatable :: ich_id(:)
      real *8, allocatable :: ts_targ(:)
      complex *16 :: zpars(3),zk
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      integer iquadtype
      complex *16, intent(out) :: wnear(nquad)
      integer, intent(out) :: iind(nquad),jind(nquad)
      real *8 rfac
      real *8, allocatable :: qwts(:)


      integer ipars
      integer ndd,ndz,ndi
      real *8 dpars,eps

      complex *16 alpha,beta,ima,val
      integer i,j,l,norder,ich
      integer ipv
      data ima/(0.0d0,1.0d0)/

      procedure (), pointer :: fker
      external h2d_slp, h2d_dlp, h2d_comb

c
c
c        initialize the appropriate kernel function
c

      iquadtype = 1
      alpha = ima
      beta = 1.0d0
      eps = 1.0d-7
      zpars(1) = zk
      zpars(2) = alpha
      zpars(3) = beta


      allocate(ich_id(ntarg),ts_targ(ntarg))

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i)
      do i=1,ntarg
        ich_id(i) = -1
        ts_targ(i) = 0.0d0
      enddo
C$OMP END PARALLEL DO     

      allocate(row_ptr(ntarg+1),col_ind(nnz))
      
      call findnear2d(cms,nch,rad_near,ndtarg,targs,ntarg,row_ptr, 
     1        col_ind)
      allocate(iquad(nnz+1)) 
      call get_iquad_rsc2d(nch,ixys,ntarg,nnz,row_ptr,col_ind,
     1         iquad)

C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,ich,norder,l)      
      do i=1,ntarg
        do j=row_ptr(i),row_ptr(i+1)-1
          ich = col_ind(j)
          norder = norders(ich)
          do l=1,norder
            iind(iquad(j)+l-1) = i
            jind(iquad(j)+l-1) = ixys(ich)+l-1 
          enddo
        enddo
      enddo
C$OMP END PARALLEL DO      

      wnear = 0
      ndz = 3
      ndi = 0
      ndd = 0
      if(iquadtype.eq.1) then
        fker => h2d_comb
        if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
          fker=>h2d_slp
        else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
          fker=>h2d_dlp
        endif
        ipv = 0
        

        call zgetnearquad_adap_guru2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,ntarg,targs,
     1     ich_id,ts_targ,eps,ipv,fker,ndd,dpars,ndz,zpars,
     1     ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,wnear)
      endif

      allocate(qwts(npts))
      call get_qwts2d(nch,norders,ixys,iptype,npts,
     1        srcvals,qwts)
     
C$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(i,j,l,val)
      do l=1,nquad
        i = iind(l)
        j = jind(l)
        call h2d_comb(srcvals(1,j),ndtarg,targs(1,i),ndd,dpars,
     1    ndz,zpars,ndi,ipars,val)
        wnear(l) = wnear(l) - val*qwts(j)

      enddo
C$OMP END PARALLEL DO


      return
      end


      subroutine getoncurvequad_helm_comb_dir_2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,
     2     eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     3     iquad,nquad,wnear)
c     
c------------------------------
c  This subroutine generates the near field quadrature 
c  for the representation 
c
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the quadrature.
c
c  Currently, the only scheme available for computing the quadrature
c  is adaptive integration
c
c  Input arguments:
c    - nch: integer
c        number of chunks
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        starting location of data on patch i
c    - iptype: integer(nch)
c        type of patch
c        iptype = 1 -> chunk discretized with Gauss Legendre nodes
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: real *8 (6,npts)
c        basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
c        if(iptype.eq.1) then basis = legendre polynomials
c    - srcvals: real *8 (8,npts)
c        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: real *8 (ndtarg,ntarg)
c        target information, the first two components must be
c        xy coordinates
c    - ich_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - ts_targ: real *8 (ntarg)
c        local t coordinate of chunk if on-surface
c    - eps: real *8
c        precision requested
c    - zpars: complex *16 (3)
c        kernel parameters (Referring to formula (1))
c        zpars(1) = k 
c        zpars(2) = alpha
c        zpars(3) = beta
c    - iquadtype - integer
c        quadrature type
c        iquadtype = 1, use adaptive quadrature for everything 
c    - nnz: integer
c        number of source patch-> target interactions in the near field
c    - row_ptr: integer(ntarg+1)
c        row_ptr(i) is the pointer to col_ind array where list of 
c        relevant source patches for target i start
c    - col_ind: integer (nnz)
c        list of source patches relevant for all targets, sorted
c        by the target number
c    - iquad: integer(nnz+1)
c        location in wnear array where quadrature for col_ind(i)
c        starts
c    - nquad: integer
c        number of entries in wnear
c
c  Output
c
c    - wnear: complex *16(nquad)
c        the desired near field quadrature
c               
c

      implicit none 
      integer, intent(in) :: nch,norders(nch),npts,nquad
      integer, intent(in) :: ixys(nch+1),iptype(nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      integer, intent(in) :: iquadtype,adjs(2,nch)
      complex *16, intent(in) :: zpars(3)
      integer, intent(in) :: nnz
      integer, intent(in) :: row_ptr(npts+1),col_ind(nnz),iquad(nnz+1)
      complex *16, intent(out) :: wnear(nquad)


      integer ipars
      integer ndd,ndz,ndi,ier,nd8
      real *8 dpars

      complex *16 alpha,beta
      integer i,j
      integer ipv

      procedure (), pointer :: fker
      external h2d_slp, h2d_dlp, h2d_comb

c
c
c        initialize the appropriate kernel function
c

      alpha = zpars(2)
      beta = zpars(3)

      ndz = 3
      ndi = 0
      ndd = 0
      if(iquadtype.eq.1) then
        fker => h2d_comb
        if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
          fker=>h2d_slp
        else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
          fker=>h2d_dlp
        endif
        ipv = 0
        

        call zgetoncurvequad_ggq2d(nch,norders,ixys,
     1       iptype,npts,srccoefs,srcvals,adjs,
     1       eps,ipv,fker,ndd,dpars,ndz,zpars,
     1       ndi,ipars,nnz,row_ptr,col_ind,iquad,nquad,wnear,ier)
      endif


      if(abs(alpha).ge.1.0d-16.and.abs(beta).lt.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*alpha
        enddo
C$OMP END PARALLEL DO        
      else if(abs(alpha).lt.1.0d-16.and.abs(beta).ge.1.0d-16) then
C$OMP PARALLEL DO DEFAULT(SHARED)        
        do i=1,nquad
          wnear(i) = wnear(i)*beta
        enddo
C$OMP END PARALLEL DO        
      endif

      return
      end
c
c
c
      
c
c
      subroutine lpcompoc_helm_comb_dir_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,adjs,
     2     eps,zpars,sigma,pot)
c
cf2py intent(in) nch,norders,ixys,iptype,npts,srccoefs,srcvals
cf2py intent(in) ndtarg,ntarg,targs,ich_id,ts_targ,eps,zpars
cf2py intent(in) sigma
cf2py intent(out) pot
c
c
c------------------------------
c  This subroutine evaluates the layer potential for the representation 
c
c
c  .. math ::
c  
c      u = (\alpha \mathcal{S}_{k} + \beta \mathcal{D}_{k})
c
c  Note: For targets on the boundary, this routine only computes
c  the principal value part, the identity term corresponding to the jump
c  in the layer potential is not included in the layer potential.
c
c
c  Input arguments:
c
c    - nch: integer
c        number of patches
c    - norders: integer(nch)
c        order of discretization on each patch 
c    - ixys: integer(nch+1)
c        ixys(i) denotes the starting location in srccoefs,
c        and srcvals array where information for patch i begins
c    - iptype: integer(nch)
c        type of patch
c    - npts: integer
c        total number of discretization points on the boundary
c    - srccoefs: double precision (6,npts)
c        koornwinder expansion coefficients of x, $\partial_{u} x$,
c        and $\partial_{v} x$. 
c    - srcvals: double precision (8,npts)
c        x, $\partial_{u} x$, $\partial_{v} x$, and $n$ sampled at
c        discretization nodes
c    - ndtarg: integer
c        leading dimension of target array
c    - ntarg: integer
c        number of targets
c    - targs: double precision (ndtarg,ntarg)
c        target information
c    - ich_id: integer(ntarg)
c        id of patch of target i, id = -1, if target is off-surface
c    - ts_targ: double precision (2,ntarg)
c        local uv coordinates on patch if on surface, otherwise
c        set to 0 by default
c    - eps: double precision
c        precision requested
c    - zpars: double complex (3)
c        kernel parameters (Referring to formula above)
c        zpars(1) = k 
c        zpars(2) = $\alpha$
c        zpars(3) = $\beta$
c     - sigma: double complex(npts)
c         density for layer potential
c
c  Output arguments
c    - pot: double complex(ntarg)
c        layer potential evaluated at the target points
c
c-----------------------------------
c
      implicit none
      integer, intent(in) :: nch,npts
      integer, intent(in) :: norders(nch),ixys(nch+1)
      integer, intent(in) :: iptype(nch),adjs(2,nch)
      real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts),eps
      complex *16, intent(in) :: zpars(3)
      complex *16, intent(in) :: sigma(npts)

      complex *16, intent(out) :: pot(npts)


      integer nptso,nnz,nquad,nd8


      integer nover,npolso
      integer norder,npols
      integer, allocatable :: row_ptr(:),col_ind(:),iquad(:)
      complex *16, allocatable :: wnear(:)

      real *8, allocatable :: srcover(:,:),wover(:)
      integer, allocatable :: ixyso(:),novers(:)

      real *8, allocatable :: cms(:,:),rads(:),rad_near(:) 

      integer i,j,jpatch,jquadstart,jstart

      integer ipars
      real *8 dpars,timeinfo(10),t1,t2,omp_get_wtime


      real *8 ttot,done,pi
      real *8 rfac,rfac0
      real *8 over4pi
      integer iptype_avg,norder_avg
      integer ikerorder, iquadtype,npts_over
      data over4pi/0.07957747154594767d0/


      iptype_avg = floor(sum(iptype)/(nch+0.0d0))
      norder_avg = floor(sum(norders)/(nch+0.0d0))

      call get_rfac2d(norder_avg,iptype_avg,rfac) 
      allocate(cms(2,nch),rads(nch),rad_near(nch))

      call get_centroid_rads2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,cms,rads)

C$OMP PARALLEL DO DEFAULT(SHARED) 
      do i=1,nch
        rad_near(i) = rads(i)*rfac
      enddo
C$OMP END PARALLEL DO      

c
c    find near quadrature correction interactions
c

      nd8=8
      call findnear2dmem(cms,nch,rad_near,nd8,srcvals,npts,nnz)

      allocate(row_ptr(npts+1),col_ind(nnz))
      
      call findnear2d(cms,nch,rad_near,nd8,srcvals,npts,row_ptr, 
     1     col_ind)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,
     1     iquad)

      ikerorder = -1
      if(abs(zpars(3)).gt.1.0d-16) ikerorder = 0

c
c    estimate oversampling for far-field, and oversample geometry
c

      allocate(novers(nch),ixyso(nch+1))
      
c     call get_far_order2d(eps,nch,norders,ixys,iptype,cms,
c     1    rads,npts,srccoefs,ndtarg,ntarg,targs,ikerorder,zpars(1),
c     2    nnz,row_ptr,col_ind,rfac,novers,ixyso)
      
      do i=1,nch
         novers(i) = norders(i)
         ixyso(i) = ixys(i)
      enddo
      ixyso(nch+1) = ixys(nch+1)
      npts_over = ixyso(nch+1)-1
      
      allocate(srcover(8,npts_over),wover(npts_over))

      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1   srccoefs,srcvals,novers,ixyso,npts_over,srcover)

      call get_qwts2d(nch,novers,ixyso,iptype,npts_over,
     1        srcover,wover)


c
c   compute near quadrature correction
c
      nquad = iquad(nnz+1)-1
      allocate(wnear(nquad))
      
C$OMP PARALLEL DO DEFAULT(SHARED)      
      do i=1,nquad
        wnear(i) = 0
      enddo
C$OMP END PARALLEL DO    


      iquadtype = 1

      call getoncurvequad_helm_comb_dir_2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,
     1     eps,zpars,iquadtype,nnz,row_ptr,col_ind,
     1     iquad,nquad,wnear)
c
c
c   compute layer potential
c
      call lpcomp_helm_comb_dir_addsub_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,nd8,npts,srcvals,
     2     eps,zpars,nnz,row_ptr,col_ind,iquad,nquad,wnear,
     3     sigma,novers,npts_over,ixyso,srcover,wover,pot)



      return
      end
c
c
c
