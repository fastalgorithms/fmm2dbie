!
!  This file contains the following user callable subroutines
!
!  get_centroid_rads2d - compute the centroid and radius of bounding
!    disk for a collection of chunks
!  
!  oversample_geom2d - oversample x,y,dxdt,dydt,dxdt2,dydt2, and normals
!    on the surface
!   
!  oversample_fun_curve2d - oversample a function defined on a curve
!
!  get_qwts2d - compute smooth quadrature weights on the curve
!  
!  get_chunk_id_ts - for all boundary points, get chunk id and local
!    t coordinate on chunk
!  
!
!

subroutine get_centroid_rads2d(nch,norders,ixys,iptype,npts,srccoefs, &
  srcvals,cms,rads)
!
!  This subroutine computes the centroid of each patch and the radius
!  of the bounding disk centered at the chunk. 
!
!  The centroid is computed as the arithmatic mean of the x,y 
!  coordinates of the discretization points on the chunks
!
!  and the bounding radius is the maximum radius of the distance to
!  the end points + all discretization points
!
!  Input arguments
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location in srccoefs,srcvals array
!        where info for each chunk begins
!    - iptype: integer(nch)
!        type of chunk
!        iptype = 1, chunk discretized with Gauss-Legendre
!        nodes and not at a corner
!    - npts: integer
!        total number of discretization points
!    - srccoefs: real *8 (6,npts)
!        Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!        if(iptype.eq.1) then basis = Legendre polynomials
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization points
!
!  Output arguments:
!    - cms: real *8 (2,nch)
!        x,y coordinates of the centroids
!    - rads: real *8 (nch)
!        radius of bounding disk centered at centroid which contains
!        the two end points of the patch
!
!
  implicit none
  integer, intent(in) :: nch,norders(nch),ixys(nch+1),iptype(nch)
  integer, intent(in) :: npts
  real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts)
  real *8, intent(out) :: cms(2,nch),rads(nch)
  
  integer i

  do i=1,nch
    if(iptype(i).eq.1) then
      call get_centroid_rads2d_lege(norders(i),srccoefs(1,ixys(i)), &
        srcvals(1,ixys(i)),cms(1,i),rads(i))
    endif
  enddo

end subroutine get_centroid_rads2d
!
!
!
!
!
subroutine get_centroid_rads2d_lege(norder,srccoefs,srcvals,cm,rad)
!
!  This subroutine computes the centroid and bounding radius
!  of a Gauss-Legendre patch in 2d
!
!  Input arguments:
!    - norder: integer
!        order of discretization of Gauss Legendre patch
!    - srccoefs: real *8 (6,norder)
!        Gauss legendre coefficients of x,y,dxdt,dydt,d2xdt2,d2ydt2
!        of the patch
!    - srcvals: real *8 (6,norder)
!        x,y,dxdt,dydt,d2xdt2,d2ydt2,rnx,rny sampled at gauss-legendre
!        nodes
! 
!  Output arguments:
!    - cm: real *8 (2)
!        xy coordinates of the centroid of the patch
!    - rad: real *8
!        radius of bounding disk which contains both end points
!
!----------------------------------------------        
  implicit none
  integer, intent(in) :: norder
  real *8, intent(in) :: srccoefs(6,norder),srcvals(8,norder)

  real *8, intent(out) :: cm(2)
  real *8, intent(out) :: rad

  real *8 ts(norder),wts(norder),umat,vmat,ra,rr,ds
  real *8 pols(norder,2),xy(2,2),t1(2)
  integer itype,i,j

  itype = 1
  call legeexps(itype,norder,ts,umat,vmat,wts)

  cm(1) = 0
  cm(2) = 0
  ra = 0
  do i=1,norder
    ds = sqrt(srcvals(3,i)**2 + srcvals(4,i)**2)*wts(i)
    ra = ra + ds
    cm(1) = cm(1) + srcvals(1,i)*ds
    cm(2) = cm(2) + srcvals(2,i)*ds
  enddo
  cm(1) = cm(1)/ra
  cm(2) = cm(2)/ra

  t1(1) = -1.0d0
  t1(2) = 1.0d0

  do i=1,2
    call legepols(t1(i),norder-1,pols(1,i))
    xy(1,i) = 0
    xy(2,i) = 0
    do j=1,norder
      xy(1,i) = xy(1,i) + srccoefs(1,j)*pols(j,i)
      xy(2,i) = xy(2,i) + srccoefs(2,j)*pols(j,i)
    enddo
  enddo
!
!  compute radius of bounding disk
!
  rad = 0
  do i=1,2
    rr = sqrt((cm(1)-xy(1,i))**2 + (cm(2)-xy(2,i))**2)
    if(rr.gt.rad) rad = rr
  enddo

end subroutine get_centroid_rads2d_lege
!
!
!
!
!


subroutine oversample_geom2d(nch,norders,ixys,iptype,npts,srccoefs, &
  srcvals,nfars,ixyso,nptso,srcover)
!
!  This subroutine oversamples geometry information
!  given expansion coeffs of geometry info on chunks
!
!
!  Extreme care must be taken when an oversampled node is identical
!  to a discretization node. The routine ensures that in this event
!  the source information is copied over from the original
!  discretization
!
!  Input arguments:
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location in srccoefs,srcvals array
!        where info for each chunk begins
!    - iptype: integer(nch)
!        type of chunk
!        iptype = 1, chunk discretized with Gauss-Legendre
!        nodes and not at a corner
!    - npts: integer
!        total number of discretization points
!    - srccoefs: real *8 (6,npts)
!        Basis coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!        if(iptype.eq.1) then basis = Legendre polynomials
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization points
!    - nfars: integer(nch)
!        oversampled order of discretization for each patch
!    - ixyso: integer(nch+1)
!        starting location of oversampled points in 
!        srcover
!    - nptso: integer
!        total number of oversampled points
!
!  Output arguments:
!    - srcover: real *8 (8,nptso)
!        oversampled geometry information
!---------------------------------

  implicit none
  integer, intent(in) :: nch,norders(nch),ixys(nch+1)
  integer, intent(in) :: iptype(nch),ixyso(nch+1)
  integer, intent(in) :: nfars(nch)
  integer, intent(in) :: npts,nptso
  real *8, intent(in) :: srccoefs(6,npts),srcvals(8,npts)
  real *8, intent(out) :: srcover(8,nptso)


  integer nfar,norder
  real *8, allocatable :: ts(:),umat,vmat,wts
  real *8, allocatable :: pmat(:,:)
  character *1 transa,transb
  real *8 alpha,beta

  integer i,istart,istarto,j,jpt,npols,l,n1,n2,ipt,i1,i2
  integer itype, nordermax, nfarmax, nfarfound, ifnfarsame
  real *8 rr

  transa = 'n'
  transb = 'n'
  alpha = 1.0d0
  beta = 0.0d0

  ifnfarsame = 1
  nfarfound = nfars(1)
  nfarmax = 0
  nordermax = 0
  do i = 1,nch
     if (nfars(i) .ne. nfarfound) ifnfarsame = 0
     nfarmax = max(nfarmax,nfars(i))
     nordermax = max(nordermax,norders(i))
  enddo

  allocate(ts(nfarmax),pmat(nordermax,nfarmax))

  if (ifnfarsame .eq. 1) then 
     itype = 0
     call legeexps(itype,nfarmax,ts,umat,vmat,wts)
     do j=1,nfarmax
        call legepols(ts(j),nordermax-1,pmat(1,j))
     enddo
  endif
     
  do i=1,nch
     nfar = nfars(i)
     norder = norders(i)
     if(nfar.ne.norder) then
        if(iptype(i).eq.1) then
           if (ifnfarsame .ne. 1) then
              itype = 0
              call legeexps(itype,nfar,ts,umat,vmat,wts)
              do j=1,nfar
                 call legepols(ts(j),norder-1,pmat(1,j))
              enddo
           endif
        endif

        istart = ixys(i)
        istarto = ixyso(i)
        call dgemm(transa,transb,6,nfar,norder,alpha, &
             srccoefs(1,istart),6,pmat,nordermax,beta,srcover(1,istarto),8)
        do j=1,nfar
           jpt = istarto + j-1
           rr = sqrt(srcover(3,jpt)**2 + srcover(4,jpt)**2)
           srcover(7,jpt) = srcover(4,jpt)/rr
           srcover(8,jpt) = -srcover(3,jpt)/rr
        enddo
        if(iptype(i).eq.1) then
           n1 = mod(norder,2)
           n2 = mod(nfar,2)
           if(n1.eq.1.and.n2.eq.1) then
              i1 = (norder+1)/2
              i2 = (nfar+1)/2
              ipt = ixys(i)+i1-1
              jpt = ixyso(i)+i2-1
              do l=1,8
                 srcover(l,jpt) = srcvals(l,ipt)
              enddo
           endif
        endif

     else
        istart = ixys(i)
        istarto = ixyso(i)
        npols = ixys(i+1)-ixys(i)
        do j=1,npols
           do l=1,8
              srcover(l,j+istarto-1) = srcvals(l,j+istart-1)
           enddo
        enddo
     endif
  enddo

end subroutine oversample_geom2d
!
!
!
!
!
subroutine oversample_fun_curv2d(nd,nch,norders,ixys,iptype,npts, &
  u,nfars,ixyso,nptso,uover)
!
!  This subroutine oversamples a collection of functions 
!  defined on a curve
!
!  Input arguments
!    - nd: integer
!        number of functions
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location in srccoefs,srcvals array
!        where info for each chunk begins
!    - iptype: integer(nch)
!        type of chunk
!        iptype = 1, chunk discretized with Gauss-Legendre
!        nodes and not at a corner
!    - npts: integer
!        total number of discretization points
!    - u: real *8 (nd,npts)
!        function values at discretization nodes
!    - nfars: integer(nch)
!        oversampled order of discretization for each patch
!    - ixyso: integer(nch+1)
!        starting location of oversampled points in 
!        srcover
!    - nptso: integer
!        total number of oversampled points
!
!  Output arguments:
!    - ufar: real *8 (nd,nptso)
!        oversampled values of functions
  implicit none
  integer, intent(in) :: nd,nch,norders(nch),ixys(nch+1)
  integer, intent(in) :: iptype(nch),npts,nptso,nfars(nch)
  integer, intent(in) :: ixyso(nch+1)
  real *8, intent(in) :: u(nd,npts)
  real *8, intent(out) :: uover(nd,nptso)

  integer i,istart,istarto,npols,npolso

  do i=1,nch
    istart = ixys(i)
    npols = ixys(i+1) - ixys(i)
    istarto = ixyso(i)
    npolso = ixyso(i+1) - ixyso(i)
    if(iptype(i).eq.1) &
      call oversample_fun_lege(nd,norders(i),u(1,istart),nfars(i),&
        uover(1,istarto))
  enddo

end subroutine oversample_fun_curv2d
!
!
!
!
!
subroutine oversample_fun_lege(nd,norder,u,nover,uover)
!
!  This subroutine oversamples a collection of functions
!  defined on [-1,1] sampled at legendre nodes
!
!  Input arguments
!    - nd: integer
!        number of functions
!    - norder: integer
!        order of discretization 
!    - u: real *8 (nd,norder)
!        function values at discretization nodes
!    - nover: integer
!        oversampled order of discretization 
!
!  Output arguments       
!    - uover: real *8 (nd,nover)
!        oversampled function values
!--------------------------        
  implicit none
  integer, intent(in) :: nd,norder,nover
  real *8, intent(in) :: u(nd,norder)
  real *8, intent(out) :: uover(nd,nover)

  real *8, allocatable :: pmat(:,:),pols(:),ucoefs(:,:)
  real *8 ts0(norder),umat0(norder,norder),vmat0(norder,norder)
  real *8 wts0(norder)
  real *8 ts(nover),umat,vmat,wts
  character *1 transa,transb
  real *8 alpha,beta
  integer itype,i

  itype = 2
  call legeexps(itype,norder,ts0,umat0,vmat0,wts0)

  itype = 0
  call legeexps(itype,nover,ts,umat,vmat,wts)
  allocate(ucoefs(nd,norder))

  transa = 'n'
  transb = 't'
  alpha = 1.0d0
  beta = 0.0d0
  call dgemm(transa,transb,nd,norder,norder,alpha,u,nd,umat0,norder, &
    beta,ucoefs,nd)

  allocate(pmat(norder,nover))
  do i=1,nover
    call legepols(ts(i),norder-1,pmat(1,i))
  enddo

  transa = 'n'
  transb = 'n'
  call dgemm(transa,transb,nd,nover,norder,alpha,ucoefs,nd,pmat,norder, &
    beta,uover,nd)


end subroutine oversample_fun_lege
!
!
!
!
!



subroutine get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,qwts)
!
!  Compute quadrature weights for integrating smooth functions
!  on the curve
!
!  Input arguments
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location in srccoefs,srcvals array
!        where info for each chunk begins
!    - iptype: integer(nch)
!        type of chunk
!        iptype = 1, chunk discretized with Gauss-Legendre
!        nodes and not at a corner
!    - npts: integer
!        total number of discretization points
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization points
!
!  Output arguments:
!    - qwts: real *8 (npts)
!        quadrature weights for integrating smooth functions
!
  implicit none
  integer, intent(in) :: nch,norders(nch),ixys(nch+1),iptype(nch),npts
  real *8, intent(in) :: srcvals(8,npts)
  real *8, intent(out) :: qwts(npts)

  real *8, allocatable :: wts(:),ts(:)
  real *8 u,v,ds
  integer itype,istart,ich,i,j,k,norder,norderfound,nordermax
  integer ifallsame

  norderfound = norders(1)
  nordermax = norders(1)
  ifallsame=1
  do i = 2,nch
     norder=norders(i)
     if (norder .ne. norderfound) ifallsame=0
     nordermax=max(nordermax,norder)
  enddo

  allocate(wts(nordermax),ts(nordermax))
  if (ifallsame .eq. 1) then
     itype = 1
     call legeexps(itype,norderfound,ts,u,v,wts)
  endif
  
  do ich=1,nch
    istart = ixys(ich)
    k = ixys(ich+1)-ixys(ich)
    if(iptype(ich).eq.1 .and. ifallsame .ne. 1) then
      itype = 1
      call legeexps(itype,k,ts,u,v,wts)
    endif
    do j=1,k
      ds = sqrt(srcvals(3,istart+j-1)**2 + srcvals(4,istart+j-1)**2)
      qwts(istart+j-1) = ds*wts(j)
   enddo
  enddo


end subroutine get_qwts2d
!
!
!
!
!
subroutine get_chunk_id_ts(nch,norders,ixys,iptype,npts,ich_id,ts_pts)
!
!  For each boundary discretization point, set chunk id
!  and local t coordinate
!
!  Input arguments
!    - nch: integer
!        number of chunks
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location in srccoefs,srcvals array
!        where info for each chunk begins
!    - iptype: integer(nch)
!        type of chunk
!        iptype = 1, chunk discretized with Gauss-Legendre
!        nodes and not at a corner
!    - npts: integer
!        total number of discretization points
!  Output arguments
!    - ich_id: integer(npts)
!        chunk number for source discretization points
!    - ts_pts: real *8 (npts)
!        local t coordinates on chunk
!      
  implicit none
  integer, intent(in) :: nch,norders(nch),ixys(nch+1),iptype(nch),npts
  integer, intent(out) :: ich_id(npts)
  real *8, intent(out) :: ts_pts(npts)

  integer ich,itype,j,k
  real *8 u,v,w
  
  do ich=1,nch
    do j=ixys(ich),ixys(ich+1)-1
      ich_id(j) = ich
    enddo
    k = ixys(ich+1)-ixys(ich)
    if(iptype(ich).eq.1) then
      itype = 0
      call legeexps(itype,k,ts_pts(ixys(ich)),u,v,w)
    endif
  enddo

end subroutine get_chunk_id_ts
!
!
subroutine get_chunk_scurv(npts,srcvals,dks)
!
!  For each boundary discretization point, compute the
! signed curvature  
!
!  Input arguments
!    - npts: integer
!        total number of discretization points
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization points
  
!  Output arguments
!    - dks: real *8(npts)
!        signed curvature at each point of geometry
!      
  implicit none
  integer, intent(in) :: npts
  real *8, intent(in) :: srcvals(8,npts)
  real *8, intent(out) :: dks(npts)

  integer i
  real *8 dx,dy,d2x,d2y,dsdt3

  do i = 1,npts
     dx = srcvals(3,i)
     dy = srcvals(4,i)     
     d2x = srcvals(5,i)
     d2y = srcvals(6,i)
     dsdt3 = sqrt(dx**2+dy**3)**3
     dks(i) = (dx*d2y-dy*d2x)/dsdt3
  enddo

end subroutine get_chunk_scurv
!
!
!
!
subroutine chunk_to_lapdlpcoef_lege(k,srcvals,targinfo,ipt,umat,rdlpcoefs)
  implicit none
  integer ipt,k
  real *8 srcvals(8,k),targinfo(2),rdlpcoefs(k),umat(k,k),fvals(k)
  real *8 x,y,rnx,rny,rself,r2,dx,dx2,dy,dy2,alpha,beta,ds
  integer i,j
 
  if(ipt.ge.1) then
    dx = srcvals(3,ipt)
    dy = srcvals(4,ipt)
    dx2 = srcvals(5,ipt)
    dy2 = srcvals(6,ipt)
    ds = sqrt(dx**2 + dy**2)
  endif

  do i=1,k
    x = targinfo(1)-srcvals(1,i)
    y = targinfo(2)-srcvals(2,i)
    rnx = srcvals(7,i)
    rny = srcvals(8,i)
    r2 = x**2 + y**2
    
    if(i.ne.ipt) then
       fvals(i) = (x*rnx + y*rny)/r2
    else
      fvals(i) = -(dx*dy2-dy*dx2)/2/ds**3 
    endif
  enddo
  alpha = 1.0d0
  beta = 0.0d0

  call dgemv('n',k,k,alpha,umat,k,fvals,1,beta,rdlpcoefs,1)
  return
end subroutine chunk_to_lapdlpcoef_lege
!
!
!
!
!
subroutine chunk_to_lapspcoef_lege(k,srcvals,targinfo,ipt,umat,rspcoefs)
  implicit none
  integer ipt,k
  real *8 srcvals(8,k),targinfo(8),rspcoefs(k),umat(k,k),fvals(k)
  real *8 x,y,rnx,rny,rself,r2,dx,d2x,dy,d2y,alpha,beta,ds
  integer i,j
  
  if(ipt.ge.1) then
    dx = srcvals(3,ipt)
    dy = srcvals(4,ipt)
    d2x = srcvals(5,ipt)
    d2y = srcvals(6,ipt)
    ds = sqrt(dx**2 + dy**2)
  endif
  rnx = targinfo(7)
  rny = targinfo(8)

  do i=1,k
    x = targinfo(1)-srcvals(1,i)
    y = targinfo(2)-srcvals(2,i)
    r2 = x**2 + y**2
    
    if(i.ne.ipt) then
       fvals(i) = (x*rnx + y*rny)/r2
    else
      fvals(i) = (dx*d2y-dy*d2x)/2/ds**3 
    endif
  enddo
  alpha = 1.0d0
  beta = 0.0d0

  call dgemv('n',k,k,alpha,umat,k,fvals,1,beta,rspcoefs,1)
  return
end subroutine chunk_to_lapspcoef_lege



subroutine chunk_to_ldlp_sp_xint(k,kg,ts,tsg,srcvals,srcvalsg,umat,xint, &
    rdlpcoefs,rspcoefs)
  implicit real *8 (a-h,o-z)
  integer k
  real *8 srcvals(8,k),xint(k,k,kg),rdlpcoefs(k,kg),rspcoefs(k,kg)
  real *8 ts(k),umat(k,k),tsg(kg),srcvalsg(8,kg)
  real *8, allocatable :: srcvals_new(:,:),dxc(:),d2xc(:),dyc(:)
  real *8, allocatable :: d2yc(:),srcvalsg_new(:)


  allocate(srcvals_new(8,k),dxc(k),d2xc(k),dyc(k),d2yc(k))
  allocate(srcvalsg_new(8))
  alpha = 1.0d0
  beta = 0.0d0
  do ipt=1,kg
    srcvals_new(1:8,1:k) = 0
    srcvalsg_new(1:8) = 0
    do i=1,k
      srcvals_new(5,i) = srcvals(5,i)*srcvalsg(7,ipt) + & 
          srcvals(6,i)*srcvalsg(8,ipt)
      srcvals_new(6,i) = -srcvals(5,i)*srcvalsg(8,ipt) + &
          srcvals(6,i)*srcvalsg(7,ipt)
    enddo
    srcvalsg_new(1) = 0 
    srcvalsg_new(2) = 0 
    srcvalsg_new(3) = 0
    srcvalsg_new(4) = sqrt(srcvalsg(3,ipt)**2 + srcvalsg(4,ipt)**2)
    srcvalsg_new(5) = srcvalsg(5,ipt)*srcvalsg(7,ipt) + &
        srcvalsg(6,ipt)*srcvalsg(8,ipt)
    srcvalsg_new(6) = -srcvalsg(6,ipt)*srcvalsg(8,ipt) + &
        srcvalsg(2,ipt)*srcvalsg(7,ipt)
    srcvalsg_new(7) = 1.0d0
    srcvalsg_new(8) = 0.0d0


      
    call dgemv('n',k,k,alpha,umat,k,srcvals_new(5,1:k),1,beta,d2xc,1)
    call dgemv('n',k,k,alpha,umat,k,srcvals_new(6,1:k),1,beta,d2yc,1)
      
    do inode=1,k
      do l=1,k
        srcvals_new(3,inode) = srcvals_new(3,inode) + & 
           xint(l,inode,ipt)*d2xc(l)
        srcvals_new(4,inode) = srcvals_new(4,inode) + &
           xint(l,inode,ipt)*d2yc(l)
      enddo
    enddo

    call dgemv('n',k,k,alpha,umat,k,srcvals_new(3,1:k),1,beta,dxc,1)
    call dgemv('n',k,k,alpha,umat,k,srcvals_new(4,1:k),1,beta,dyc,1)
      
    do inode=1,k
      do l=1,k
        srcvals_new(1,inode) = srcvals_new(1,inode) + &
           xint(l,inode,ipt)*dxc(l)
        srcvals_new(2,inode) = srcvals_new(2,inode) + &
           xint(l,inode,ipt)*dyc(l)  
      enddo
    enddo

    iptuse = -1
    do inode=1,k
      srcvals_new(4,inode) = srcvals_new(4,inode)+ &
         sqrt(srcvalsg(3,ipt)**2 + srcvalsg(4,ipt)**2)
      srcvals_new(2,inode) = srcvals_new(2,inode) + &
        sqrt(srcvalsg(3,ipt)**2 + srcvalsg(4,ipt)**2)*(ts(inode)-tsg(ipt))
      if(abs(ts(inode)-tsg(ipt)).le.1.0d-16) iptuse = inode
      ds = sqrt(srcvals_new(3,inode)**2 + srcvals_new(4,inode)**2)
      srcvals_new(7,inode) = srcvals_new(4,inode)/ds
      srcvals_new(8,inode) = -srcvals_new(3,inode)/ds
    enddo
    call chunk_to_lapdlpcoef_lege(k,srcvals_new,srcvalsg_new, &
       iptuse,umat,rdlpcoefs(1,ipt))
    call chunk_to_lapspcoef_lege(k,srcvals_new,srcvalsg_new, &
       iptuse,umat,rspcoefs(1,ipt))
  enddo

  return
end subroutine
