subroutine chunk_interior(nch,norders,ixys,iptype,npts, &
     srcvals,srccoefs,targs,ndt,nt,inflag)
  !
  ! decide whether or not a point is inside a chunk
  ! geometry. The curve should be positively oriented
  ! (traversed clockwise so that the interior is on the
  ! left). This routine should be reasonably efficient
  ! It has only been tested on single closed curves.    
  !
  ! Input:
  !
  ! CHUNK GEOMETRY FORMAT
  !
  ! nch, norders, ixys, iptype, npts, srcvals, srccoefs
  !
  ! OTHER INPUT
  !
  ! targs - real *8 (ndt,nt) array, targs(1:2,i) is the
  !          location of the ith point to be identified
  !          as inside or outside the domain
  !
  ! ndt - integer, leading dimension of targs array
  !
  ! Output:
  !
  ! inflag - integer(nt), inflag(i) = 1 if targs(1:2,i)
  !               is inside the domain, inflag(i) = -1
  !               otherwise
  !
  implicit real *8 (a-h,o-z)
  real *8 :: srcvals(8,npts), srccoefs(6,npts), targs(ndt,nt)
  integer :: inflag(nt), ixys(nch+1), iptype(nch), norders(nch)
  complex *16, allocatable :: dipstr(:), pot(:), pottarg(:)
  real *8, allocatable :: wts(:),src(:,:),srcwid(:), &
       cms(:,:),rads(:),targ1(:,:),ts_targ(:),dist_targ(:), &
       cf(:,:)
  integer, allocatable :: ich_id(:), icheck(:)

  real *8 timeinfo(8)
  complex *16 :: eye, twopieye
  data eye /(0.0d0,1.0d0)/

  pi = 4*atan(1.0d0)
  
  twopieye = 2*pi*eye


  ! first pass is Cauchy's ID
  
  allocate(wts(npts))
  call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,wts)

  allocate(dipstr(npts),src(2,npts))

  do i = 1,npts
     src(1,i) = srcvals(1,i)
     src(2,i) = srcvals(2,i)
     dipstr(i) = -wts(i)*(srcvals(3,i)+eye*srcvals(4,i))
     dipstr(i) = dipstr(i)/twopieye
     dipstr(i) = dipstr(i)/sqrt(srcvals(3,i)**2+srcvals(4,i)**2)
  enddo

  allocate(pot(npts),targ1(2,nt),pottarg(nt))

  do i = 1,nt
     targ1(1,i) = targs(1,i)
     targ1(2,i) = targs(2,i)
  enddo

  eps = 1d-12
  call cpu_time(t1)
  !$ t1 = omp_get_wtime()  
  call cfmm2d_st_d_p(eps,npts,src, &
       dipstr,pot,nt,targ1,pottarg)
  call cpu_time(t2)
  !$ t2 = omp_get_wtime()    

  ! anything within tol of 1 is inside
  ! anything within tol of 0 is outside
  tol = 1d-5

  allocate(icheck(nt))
  
  nin = 0
  nout = 0
  nchk = 0
  do i = 1,nt
     inflag(i) = 0
     if (abs(pottarg(i)-1) .lt. tol) inflag(i)=1
     if (abs(pottarg(i)) .lt. tol) inflag(i)=-1

     if (inflag(i).eq. 1) nin=nin+1
     if (inflag(i) .eq. -1) nout=nout+1
     if (inflag(i) .eq. 0) then
        nchk = nchk+1
        icheck(nchk) = i
        targ1(1,nchk) = targs(1,i)
        targ1(2,nchk) = targs(2,i)        
     endif
        
  enddo

  ! for vague points, find nearest boundary
  ! point and use normal direction test

  allocate(cms(2,nch),rads(nch),srcwid(npts))
  call get_centroid_rads2d(nch,norders,ixys,iptype,npts, &
       srccoefs,srcvals,cms,rads)

  do i = 1,nch
     do j = ixys(i),(ixys(i+1)-1)
        srcwid(j) = 2*rads(i)
     enddo
  enddo

  allocate(ich_id(nchk),ts_targ(nchk),dist_targ(nchk))

  nd2=2
  call findnearchunktarg_id_ts(nch,norders,ixys,iptype,npts, &
       srccoefs,srcvals,srcwid,nd2,nchk,targ1,ich_id,ts_targ, &
       dist_targ,timeinfo,ier)

  ! normal test

  kmax = 1
  do i = 1,nch
     kmax = max(kmax,norders(i))
  enddo

  allocate(cf(kmax,4))
  
  do i = 1,nchk
     irel = icheck(i)
     ich = ich_id(i)
     t = ts_targ(i)
     xt = targ1(1,i)
     yt = targ1(2,i)

     if (ich .le. 0) then
        ! failed, revert to best guess
        inout = 1
        if (abs(pottarg(irel)) .lt. abs(pottarg(irel)-1)) &
             inout = -1
     else
        k = norders(ich)
        kord = k-1
        jstart = ixys(ich)
        do j = 1,k
           cf(j,1)=srccoefs(1,jstart+j-1)
           cf(j,2)=srccoefs(2,jstart+j-1)
           cf(j,3)=srccoefs(3,jstart+j-1)
           cf(j,4)=srccoefs(4,jstart+j-1)           
        enddo
        call legeexev(t,x,cf(1,1),kord)
        call legeexev(t,y,cf(1,2),kord)
        call legeexev(t,dx,cf(1,3),kord)
        call legeexev(t,dy,cf(1,4),kord)
        
        c = (x-xt)*dy-(y-yt)*dx

        inout = -1
        if (c .ge. 0) inout = 1
        
     endif

     inflag(irel) = inout
  enddo
   
  
  return
end subroutine chunk_interior


