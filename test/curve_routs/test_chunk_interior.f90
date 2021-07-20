program test001
  implicit real *8 (a-h,o-z)

  call prini(6,13)
  iseed = 8675309
  d = hkrand(iseed)
  
  isuccess=0
  call testit(isuccess)

  if (isuccess .eq. 1) write(*,*) 'success!'
  
  stop
end program test001


subroutine testit(isuccess)
  implicit real *8 (a-h,o-z)
  !
  real *8, allocatable :: srcvals(:,:),srccoefs(:,:),ab(:,:)
  integer, allocatable :: norders(:),ixys(:),iptype(:),adjs(:,:)
  !
  real *8 dpars(2) 
  integer ipars(1)
  complex *16 zpars
  external fstarn_simple
  !
  real *8, allocatable :: targs(:,:)
  integer, allocatable :: inflag(:), intrue(:)
  !
  
  isuccess = 0

  done = 1
  pi = atan(done)*4

  ! discretize curve

  nchmax = 200000
  k = 14
  npts_max = nchmax*k
  allocate(srcvals(8,npts_max),srccoefs(6,npts_max),ab(2,nchmax))
  allocate(norders(nchmax),ixys(nchmax+1),iptype(nchmax))
  allocate(adjs(2,nchmax))

  irefinel = 0
  irefiner = 0
  rlmax = 1.0d10
  ta = -pi
  tb = pi
  ifclosed = 1
  rlmaxe = 0

  ndd = 2
  ndi = 1
  ndz = 0

  ! parameters sent to routine defining the curve (fstarn_simple)
  dpars(1) = 1.0d0
  dpars(2) = 0.25d0
  ipars(1) = 5

  ! discretization paratmers
  nover = 0
  nch = 0
  ier = 0
  epsc = 1.0d-6
  call chunkfunc_guru(epsc,rlmax,ifclosed,irefinel,irefiner, &
       rlmaxe,ta,tb,fstarn_simple,ndd,dpars,ndz,zpars,ndi, &
       ipars,nover,k,nchmax,nch,norders,ixys,iptype,npts, &
       srcvals,srccoefs,ab,adjs,ier)


  call prinf('ier=*',ier,1)
  if(ier.ne.0) then
     call prinf('failed to execute*',ier,1)
     stop
  endif
  call prinf('nch=*',nch,1)


  ! set up test points and call routine

  ndt=4
  nt = 100000

  allocate(inflag(nt),intrue(nt),targs(ndt,nt))

  do i = 1,nt
     t=-pi + 2*pi*hkrand(0)
     call fstarn_simple(t,ndd,dpars,ndz,zpars,ndi,ipars, &
          x,y,dxdy,dydt,dxdt2,dydt2)
     scal = 2*hkrand(0)
     intrue(i)=-1
     if (scal .le. 1) intrue(i)=1
     x=x*scal
     y=y*scal
     targs(1,i)=x
     targs(2,i)=y
  enddo
  
  call cpu_time(t1)
  call chunk_interior(nch,norders,ixys,iptype,npts, &
       srcvals,srccoefs,targs,ndt,nt,inflag)
  call cpu_time(t2)
  call prin2('time taken in near target generation=*',t2-t1,1)

  ! check values

  isuccess=1
  do i = 1,nt
     if (inflag(i).ne.intrue(i)) isuccess=0
  enddo

  return
end subroutine testit


subroutine fstarn_simple(t,ndd,dpars,ndz,zpars,ndi,ipars,x, &
     y,dxdt,dydt,dxdt2,dydt2)
!
!   curve parametrization for simple star shaped domain
!
!       r(t) = r1 + r2*cos(n*t)
!       x(t) = r(t) cos(t)
!       y(t) = r(t) sin(t)
!   
!   where r1 = dpars(1), r2 = dpars(2) and n = ipars(1)
!     
!
  implicit real *8 (a-h,o-z)
  real *8 dpars(2)
  integer ipars
  complex *16 zpars,z,zd,zd2,zexp,ima,ze1,ze2,ze3
  data ima/(0.0d0,1.0d0)/

  rsc = 1.0d0

  n = ipars
  r1 = dpars(1)*rsc
  r2 = dpars(2)*rsc

  r = r1 + r2*cos(n*t)


  drdt = -r2*(n+0.0d0)*sin(n*t)
  d2rdt2 = -r2*(n+0.0d0)*(n+0.0d0)*cos(n*t)

  x = r*cos(t)
  y = r*sin(t)

  dxdt = drdt*cos(t) - r*sin(t)
  dydt = drdt*sin(t) + r*cos(t)

  dxdt2 = d2rdt2*cos(t) - 2*drdt*sin(t) - r*cos(t) 
  dydt2 = d2rdt2*sin(t) + 2*drdt*cos(t) - r*sin(t) 

  z = (r1 + r2*cos((n+0.0d0)*t))*exp(ima*t)
  zd = -(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t) + ima*z
  zd2 = ima*zd - (n+0.0d0)**2*r2*cos((n+0.0d0)*t)*exp(ima*t) - &
       ima*(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t)
  ze1 = exp(ima*t)
  ze2 = exp(ima*t*(n+1))
  ze3 = exp(ima*t*(-n+1))
  z = r1*ze1 + 0.5d0*r2*(ze2+ze3)
  zd = ima*(r1*ze1 + 0.5d0*r2*(n+1)*ze2 + 0.5d0*r2*(-n+1)*ze3)
  zd2 = -(r1*ze1 + 0.5d0*r2*(n+1)**2*ze2 + 0.5d0*r2*(-n+1)**2*ze3)

  x = real(z)
  y = imag(z)
  dxdt = real(zd)
  dydt = imag(zd)
  dxdt2 = real(zd2)
  dydt2 = imag(zd2)


  return
end subroutine fstarn_simple
