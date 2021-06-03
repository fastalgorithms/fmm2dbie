      subroutine dchunkints_lege(eps,norder,srccoefs,ndtarg,
     1   ntarg,xytarg,nporder,fker,ndd,dpars,ndz,zpars,ndi,
     2   ipars,nqorder,nchmax,dintvals)
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
c   until an absolute tolerance of $\eps$ is reached.
c
c   Notes:
c   Unlike the analogous 3d routines, there are a few differences
c   for the 2d routines. 
c   1) we only support adaptive integration
c   as the method of computing the integrals, there aren't any
c   considerations based on distance to the targets
c   2) Proxy target locations have now been eliminated and absorbed in
c   the xytargs array
c   3) The routines are currently not optimized to reuse the
c   geometry info at multiple targets or the legendre polynomials
c   at the multiple targets
c   4) Legendre polynomial evaluation also isn't optimized and can
c   be improved by precomputing the recursion relation coefficients
c
c   Input arguments:
c     - eps: real *8
c         requested tolerance
c     - norder: integer
c         order of discretization of input chunk
c     - srccoefs: real *8 (6,norder)
c         Legendre expansion coefficients of x,y,dxdt,dydt,d2xdt2,d2ydt2
c     - ndtarg: integer
c         leading order dimension of target location array
c     - ntarg: integer
c         number of targets
c     - xytarg: real *8 (ndtarg,ntarg)
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
c     - dintvals: real *8 (nporder,ntarg)
c         computed integrals
c   
c-----------------------------------
      implicit none
      real *8, intent(in) :: eps
      integer, intent(in) :: norder,ndtarg,ntarg,nporder,ndd,ndz,ndi
      integer, intent(in) :: ipars(ndi),nqorder,nchmax
      real *8, intent(in) :: srccoefs(6,norder),dpars(ndd)
      real *8, intent(in) :: xytarg(ndtarg,ntarg)
      complex *16, intent(in) :: zpars(ndz)

      real *8, intent(out) :: dintvals(nporder,ntarg)

      real *8, allocatable :: stack(:,:)
      real *8, allocatable :: vals(:,:)
      real *8, allocatable :: value2(:),value3(:)
      real *8, allocatable :: ts0(:),w0(:)
      real *8 umat0,vmat0
      integer itype

      real *8 a,b
      integer maxrec,numit,maxdepth,nnmax,i,j,l,ier

      external fker

      itype = 1
      allocate(ts0(nqorder),w0(nqorder))
      call legeexps(itype,nqorder,ts0,umat0,vmat0,w0)

      maxdepth = 200
      allocate(stack(2,maxdepth),vals(nporder,maxdepth))
      allocate(value2(nporder),value3(nporder))

      a = -1.0d0
      b = 1.0d0

      do i=1,ntarg
        nnmax = nchmax
        do j=1,maxdepth
          stack(1,j) = 0
          stack(2,j) = 0
          do l=1,nporder
            vals(l,j) = 0
          enddo
        enddo
        
        ier = 0
        maxrec = 0
        numit = 0
        call dadinrecm(ier,stack,a,b,norder,srccoefs,ndtarg,xytarg(1,i),
     1     nporder,fker,ndd,dpars,ndz,zpars,ndi,ipars,nqorder,ts0,w0,
     2     vals,nnmax,eps,dintvals(1,i),maxdepth,maxrec,numit,value2,
     3     value3)
        print *, i,numit
      enddo

      return
      end
c
c
c
c
c
      subroutine dadinrecm(ier,stack,a,b,k,srccoefs,ndt,targ,np,fker,
     1   ndd,dpars,ndz,zpars,ndi,ipars,m,ts0,w0,vals,nnmax,
     2   eps,rints,maxdepth,maxrec,numit,value2,value3)
c
c 
c  This subroutine computes the integrals 
c      
c     \int_{a}^{b} K(x,\rho(t)) P_{n}(t) dsdt(t) dt \, ,
c  
c   where $P_{n}(t)$ are the legendre polynomials,
c   $\rho(t)$ is the parametrization of a chunk
c   dsdt is the speed of parametrization of $\rho(t)$,
c   for one target location $x$
c
c   The map $\rho(t)$ is stored as coefficients of a legendre expansion
c   $xy(t)$ along with expansions of the first and second derivative
c   information with respect to t.
c
c   The integrals are computed using adaptive integration 
c
c   Input arguments:
c     - stack: real *8 (2,maxdepth)
c         will be storing the stack of subdivisions of the interval
c         used during the adaptive integration process.
c         This is just one of the temporary work arrays
c     - a: real *8
c         left end point of the interval
c     - b: real *8
c         right end point of the interval
c     - k: integer
c         order of discretization of input chunk
c     - srccoefs: real *8 (6,k)
c         Legendre expansion coefficients of x,y,dxdt,dydt,d2xdt2,d2ydt2
c     - ndtarg: integer
c         leading order dimension of target location array
c     - targ: real *8 (ndtarg)
c         target info. The first two components xytarg(1:2,i) must be 
c         the xy coordinates of the targets
c     - np: integer
c         order of legendre polynomials to be integrated
c     - fker: function handle
c         function handle for evaluating the kernel K
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
c     - nnmax: integer
c         maximum number of chunks in adaptive integration, 
c         will be increased by factors of 2 after max is reached
c         for current chunk (Feature under construction)
c     - eps: real *8
c         tolearance requested
c     - maxdepth: integer
c         maximum depth in adaptive integration
c     - vals: real *8 (np,maxdepth)
c         work array 1
c     - value2: real *8 (np)
c         work array 2
c     - value3: real *8 (np)
c         work array 3
c    
c   Output arguments:
c     - ier: integer
c         error code.
c         * ier = 0, successful execution
c         * ier = 8, one subinterval in the subdivision
c           was smaller than (b-a)/2**maxdepth
c         * ier = 16, total number of subintervals 
c           in adaptive integration turned out to be
c           greater than nnmax
c     - rints: real *8 (np)
c         computed integrals
c     - maxrec: integer
c         maximum depth to which the recursion went at its deepest
c         point
c     - numint: integer
c         total number of intervals in the subdivsion
c     
c   
c
c
     
      implicit real *8 (a-h,o-z)
      real *8 stack(2,maxdepth),a,b
      integer k,np
      real *8 srccoefs(6,k),targ(ndt),dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 ts0(m),w0(m)
      real *8 vals(np,maxdepth),rints(np),value2(np),value3(np)
      external fker
      
      
      ier = 0
      stack(1,1) = a
      stack(2,1) = b

      call doneintmu(a,b,k,srccoefs,ndt,targ,np,fker,ndd,dpars,ndz,
     1   zpars,ndi,ipars,m,ts0,w0,vals(1,1))
      
c 
c       recursively integrate the thing
c 
      j=1
      do i=1,np
        rints(i)=0
      enddo
c 
      ier=0
      maxrec=0
      do 3000 i=1,nnmax
        numit=i
        if(j .gt. maxrec) maxrec=j
c 
c       subdivide the current subinterval
c 
         c=(stack(1,j)+stack(2,j))/2
         call doneintmu(stack(1,j),c,k,srccoefs,ndt,targ,np,fker,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value2)
         call doneintmu(c,stack(2,j),k,srccoefs,ndt,targ,np,fker,ndd,
     1      dpars,ndz,zpars,ndi,ipars,m,ts0,w0,value3)
c 
          dd=0
          do jj=1,np
            ddd=abs(value2(jj)+value3(jj)-vals(jj,j) )
            if(ddd .gt. dd) dd=ddd
          enddo
c 
          ifdone=0
          if(dd .le. eps) ifdone=1
c 
c       if the function on this subinterval has been
c       integrated with sufficient accuracy - add the
c       value to that of the global integral and move up
c       in the stack
c 
          if(ifdone  .eq. 0) goto 2000
c 
          do jj=1,np
            rints(jj)=rints(jj)+value2(jj)+value3(jj)
          enddo
          j=j-1
c 
c        if the whole thing has been integrated - return
c 
          if(j .eq. 0) return
          goto 3000
 2000 continue
c 
c       if the function on this subinterval has not been
c       integrated with sufficient accuracy - move
c       down the stack
c 
        stack(1,j+1)=stack(1,j)
        stack(2,j+1)=(stack(1,j)+stack(2,j))/2
        do jj=1,np
          vals(jj,j+1) = value2(jj)
          vals(jj,j) = value3(jj)
        enddo
c 
        stack(1,j)=(stack(1,j)+stack(2,j))/2
c 
        j=j+1
c 
c       if the depth of the recursion has become excessive - bomb
c 
        if(j .lt. maxdepth) goto 3000
        ier=8
        return
 3000 continue
      ier=16

      return
      end
c
c
c
c
c
      subroutine doneintmu(a,b,k,srccoefs,ndt,targ,np,fker,ndd,dpars,
     1   ndz,zpars,ndi,ipars,m,t,w,vals)
      implicit real *8 (a-h,o-z)
      real *8 a,b
      real *8 srccoefs(6,k),targ(ndt)
      real *8 srctmp(8),pols(k+np)
      real *8 dpars(ndd)
      complex *16 zpars(ndz)
      integer ipars(ndi)
      real *8 t(m),w(m)
      real *8 vals(np),val
      external fker
      
      do i=1,np
        vals(i) = 0
      enddo

      u = (b-a)/2
      v = (b+a)/2
      alpha = 1.0d0
      beta = 0.0d0
      nuse = max(k,np)

      do i=1,m
        tt = u*t(i)+v
        call legepols(tt,nuse-1,pols) 
        call dgemv('n',6,k,alpha,srccoefs,6,pols,1,beta,srctmp,1)
        dst = sqrt(srctmp(3)**2 + srctmp(4)**2)
        srctmp(7) = srctmp(4)/dst
        srctmp(8) = -srctmp(3)/dst
        val = 0
        call fker(srctmp,ndt,targ,ndd,dpars,ndz,zpars,ndi,
     1      ipars,val)
        do j=1,np
          vals(j) = vals(j) + dst*w(i)*val*pols(j) 
        enddo
      enddo

      do i=1,np
        vals(i) = vals(i)*u
      enddo

      return
      end

