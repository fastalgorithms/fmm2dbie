!     Copyright (C) 2019: Michael O'Neil
!     Contact: oneil@cims.nyu.edu
!      
!c This software is being released under a modified FreeBSD license
!c (see LICENSE in home directory). 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!       $Date$
!       $Revision$
!
!  The code below contains several user-callable routines
!  for the construction of refind polygons chunked up smooth
!  curves. The most useful routines are:
!
!  chunkfunc_guru - passing in a subroutine which describes a 
!             curve analytically, returns an adpative chunked
!             version given by Legendre nodes.
!             Also returns the intervals used in parameter
!             space and the adjacency information
!
!  chunklength - receives the same subroutine that chunkfunc
!             does, calculates the arclength between
!             two specified points (in parameterization space)
!
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
      subroutine chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner, &
         rlmaxe,ta,tb,funcurve,ndd,dpars,ndz,zpars,ndi,ipars,nover, &
         k,nchmax,nch,norders,ixys,iptype,npts,srcvals, &
         srccoefs,ab,adjs,ier)
!
!  Using a user-defined subroutine funcurve, split up the
!  curve into chunks. the calling sequence of funcurve
!  should be
!
!  funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars,x,y,dxdt,dydt,dxdt2,dydt2)
!  
!
!  The routine refines the chunks if 
!  1) their length is greater than  rlmax  OR 
!  2) if chunk is not resolved OR 
!  3) if chunk at an end point ta or tb is greater than rlmaxe as long as
!  irefinel or irefiner = 1
!
!  The first condition is meant to allow for sampling of the curve
!  at fixed number of points per wavelength
!
!  The second condition is to ensure that the curve is resolved
!  to desired tolerance
!
!  The third condition is to allow dyadic refinement of the discretization
!  at the end points which might be necessary for regions with corners
!  or open arcs
!  The routine will declare a chunk to be resolved
!  if the relative l2 error on the chunk for the g = curvature^2*dsdt
!  satisfies |g - ginterp|*(rlpan)/(1+|g-ginterp|*sqrt(rlpan)) \leq 
!  eps_test
!  and the relative l2 error for dsdt on the chunk is \leq eps_test
!
!  here eps_test is determined based on the user described tolerance
!  eps and the conditioning of the problem. It is given by
!  \eps_test = max(eps,2^(-48)*(1+|kmax|)(1+|kmin|)
!  where kmax and kmin are the maximum and minimum absolute values
!  of the curvature
!
!  and (rlpan) above is the length of the panel
!
!  NOTE: the routine returns chunks that are at most a factor
!   of two different in archlength than adjacent chunks
!
!  Input arguments:  
!    - eps: double precision
!        absolute precision to resolve the curve
!    - rlmax: double precision
!        maximum permissible chunk length, if chunk length >rlmax,
!        it will be refined
!    - ifclosed: integer
!        flag for indicating if curve is closed. This parameter
!        is needed to enforce 2:1 panel length criterion and
!        deciding how to compute neighbors of panels at
!        the end points
!    - irefinel, irefiner: integeer
!        flags for refining discretization at either end point
!        end point left (ta) or right (tb) will be refined
!        if irefinel = 1 or irefiner = 1
!    - rlmaxe: double precision
!        relevant only if irefinel or irefiner = 1, chunk
!        at end point will be refined if chunk length > rlmaxe
!    - ta,tb: double precision
!        assume that curve is parameterized by t \in [ta,tb)
!    - funcurve: function handle
!        function which returns x,y,dxdt,dydt,dxdt2,dydt2
!        for a given point in parameter space
!    - ndd: integer
!        length of double precision parameter array for funcurve
!    - dpars: double precision(ndd)
!        list of double precision parameters for funcurve
!    - ndz: integer
!        length of complex parameter array for funcurve
!    - dpars: double comlex(ndd)
!        list of complex parameters for funcurve
!    - ndi: integer
!        length of integer parameter array for funcurve
!    - ipars: integer(ndi)
!        list of integer parameters for funcurve
!    - nover: integer
!        post-process oversampling factor, nover <= 1 will result
!        in no change, nover=2 will split each chunk in half (with
!        respect to arclength!!)
!    - k: integer
!        number of legendre nodes per chunk
!    - nchmax: integer
!        max number of chunks in the discretization. The subroutine
!        will return with an error if the final number of chunks
!        required is more than nchmax
!
!  Note on input, all the output arrays dimensioned by nch
!  should be of size nchmax, but only the first nch (nch+1 for ixys)
!  entries will be populated. Similarly the points and coefs
!  array must be of size nchmax*k, but only the first k*nch
!  entries will be populated 
!
!  Output arguments:
!    - nch: integer
!        total number of chunks created
!    - norders: integer(nch)
!        order of discretization of each chunk
!    - ixys: integer(nch+1)
!        starting location of points on chunk i in
!        srcvals,srccoefs array
!    - iptype: integer(nch)
!        type of chunks in discretization. This discretization
!        only returns chunks of type=1
!    - npts: integer
!        total number of points in the discretization. npts = k*nch.
!    - srcvals: real *8 (8,npts)
!        x,y,dxdt,dydt,dxdt2,dydt2,rnx,rny at the discretization points
!    - srccoefs: real *8 (6,npts)
!        Legendre expansion coefficients of x,y,dxdt,dydt,dxdt2,dydt2
!    - ab: real *8 (2,nch) 
!        left and right end points in parameter space of each chunk
!    - adjs: integer(2,nch)
!        adjacency information for each chunk
!    - ier: integer
!        error code. 
!        ier = 0 => successful generation of curve
!        ier = 4 => Not sufficient memory. nchmax too small, rerun with 
!           larger nchmax
!
!
!
      implicit real *8 (a-h,o-z)
      
      ! List of calling sequence variables
      real *8, intent(in) :: eps,rlmax,rlmaxe,ta,tb,dpars(ndd)
      integer, intent(in) :: irefinel,irefiner,ndd,ndi,ndz,ipars(ndi)
      integer, intent(in) :: k,nchmax,ifclosed
      complex *16, intent(in) :: zpars(ndz)
      integer, intent(out) :: nch,norders(nchmax),ixys(nchmax), &
        iptype(nchmax),npts,ier,adjs(2,nchmax)
      real *8, intent(out) :: srcvals(8,*),srccoefs(6,*),ab(2,nchmax)

      

      external funcurve
!
!     Temporary variables
!  
      integer, allocatable :: ifprocess(:)
      real *8, allocatable :: xs(:),ws(:),u(:,:),v(:,:)
      real *8, allocatable :: xs2(:),ws2(:),u2(:,:),v2(:,:),ttest(:)
      real *8, allocatable :: fvals(:,:),finterp(:,:),ximat(:,:)
      real *8, allocatable :: finterpex(:,:),work(:),curv(:)
      real *8 rints(4),rints0(4),rinttmp(4),errs(4),errs0(4),epsuse
      real *8 rintssq(4)
      real *8 cxy(2),cxy0(2),cxytmp(2)
      complex *16 z,zd,zd2,ima
      data ima/(0.0d0,1.0d0)/
      
      done=1
      pi=4*atan(done)
      ier = 0
!
      allocate(ifprocess(nchmax))
      do i=1,nchmax
        ifprocess(i)=0
      enddo

      allocate(xs(k),ws(k),u(k,k),v(k,k))

!
!       construct legendre nodes and weights, k and 2k of them, as well
!       as the interpolation/coefficients matrices
!
      itype=2
      call legeexps(itype,k,xs,u,v,ws)

      allocate(ttest(k))
      do i=1,k
        ttest(i) = -1.0d0 + 2.0d0*(i-1.0d0)/(k-1.0d0)
      enddo

  
      

!
!       . . . start chunking
!
      ab(1,1)=ta
      ab(2,1)=tb
      nch=1
      ifdone=1
      adjs(1,1)=-1
      adjs(2,1)=-1
      nchnew=nch

      allocate(fvals(2,k))
      allocate(finterpex(2,k),finterp(2,k),ximat(k,k))
      lw = 2*k**2 + k+100
      allocate(work(lw))
      call lematrin(k,k,ttest,ximat,xs,work)
      call interpmat_1d_hels(k,ttest,xs,ximat)


      
      
      
!
!  Estimate arc length of curve
!
!
      rlcurve = 0
      call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
         ta,tb,xs,ws,rlcurve)
      rlcurve0 = rlcurve
      allocate(curv(k))
      call chunkcurv(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars,&
        ta,tb,xs,ws,curv)
      

      alpha = 1.0d0
      beta = 0.0d0

 1311 format(2x,i1,2x,i4,5(2x,e11.5))
      maxiter=52
      do ijk = 1,maxiter
        rlcurve = rlcurve0
        epsuse = eps
!        read *, i
!
!     loop through all existing chunks, if resolved store, if not split
!
        ifdone=1
        do ich=1,nchnew
!
          if (ifprocess(ich) .eq. 1) goto 4600
          ifprocess(ich)=1
!
          a=ab(1,ich)
          b=ab(2,ich)
          iprint = 0

!          if(ijk.eq.11.and.ich.eq.368) then
!            iprint = 1
!            print *, ich,a,b
!            a = a - 2*pi
!            b = b - 2*pi
!          endif

!
          rlpan = 0
          do i=1,k
            t=a+(b-a)*(1+xs(i))/2
            call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
            x,y,dx,dy,dx2,dy2)

            fvals(1,i) = sqrt(dx**2 + dy**2)
            zd = dx + ima*dy
            zd2 = dx2 + ima*dy2
            fvals(2,i) = dimag(zd2*dconjg(zd))**2/abs(zd)**5 
            rlpan = rlpan + fvals(1,i)*ws(i)*(b-a)/2
            
          enddo

          do i=1,k
            t = a+(b-a)*(1+ttest(i))/2
            call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
              x,y,dx,dy,dx2,dy2)

            finterpex(1,i) = sqrt(dx**2 + dy**2)
            zd = dx + ima*dy
            zd2 = dx2 + ima*dy2
            finterpex(2,i) = dimag(zd2*dconjg(zd))**2/abs(zd)**5
            rk = abs(imag(zd2/zd)/finterpex(1,i))+1
            if(rk.gt.rkmax) rkmax = rk
            if(rk.lt.rkmin) rkmin = rk
            
          enddo

          call dgemm('n','t',2,k,k,alpha,fvals,2,ximat,k, &
              beta,finterp,2)

          err1a = 0
          err2a = 0
          err1r = 0
          err2r = 0
          r1 = 0
          r2 = 0
          do i=1,k
            err1 = abs(finterp(1,i)-finterpex(1,i))
            if(err1.gt.err1a) err1a = err1
            err2 = abs(finterp(2,i)-finterpex(2,i))
            if(err2.gt.err2a) err2a = err2
            r1 = r1 + abs(finterpex(1,i))
            r2 = r2 + abs(finterpex(2,i))
          enddo
          err1r = err1a/r1
          err2r = err2a/r2

          err1 = min(err1a,err1r)
          err2 = min(err2a,err2r)

!            
!
!
!       . . . mark as processed and resolved if less than eps
!
!       Note we are checking both resolution and panel length
!       simultaneously since the length might be inaccurate
!       if the panel is not resolved anyway. 
!
!       There could some potential savings by doing the
!       panel check earlier which could be explored. The current
!       extra cost is just a safety measure for robustness
!
!
!
          rmsemax = 0.0d0

          if(err1.gt.epsuse) then
            goto 2800
          endif

          if(err2.gt.epsuse) then
            goto 2800
          endif


          if (rlpan.gt.rlmax) goto 2800
          goto 4600
!
 2800 continue
!
!       . . . if here, not resolved
!       divide - first update the adjacency list
!
            
          ifprocess(ich)=0
          ifdone=0

          if ((nch .eq. 1) .and. (ifclosed .gt. 0)) then
            adjs(1,nch)=2
            adjs(2,nch)=2
            adjs(1,nch+1)=1
            adjs(2,nch+1)=1
          endif
!
          if ((nch .eq. 1) .and. (ifclosed .le. 0)) then
            adjs(1,nch)=-1
            adjs(2,nch)=2
            adjs(1,nch+1)=1
            adjs(2,nch+1)=-1
          endif
!
          if (nch .gt. 1) then
            iold2=adjs(2,ich)
            adjs(2,ich)=nch+1
            if (iold2 .gt. 0) adjs(1,iold2)=nch+1
            adjs(1,nch+1)=ich
            adjs(2,nch+1)=iold2
          endif

!
!       now update the endpoints in ab
!
          ab(1,ich)=a
          ab(2,ich)=(a+b)/2
!
          nch=nch+1
          if (nch .gt. nchmax) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            ier = 4
            return
          endif
!
          ab(1,nch)=(a+b)/2
          ab(2,nch)=b
!
!         Update rlcurve0
!
          rlcurve0 = rlcurve0 - rlpan
          a = ab(1,ich)
          b = ab(2,ich)
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,b,xs,ws,rlpan)
          rlcurve0 = rlcurve0 + rlpan

          a = ab(1,nch)
          b = ab(2,nch)
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,b,xs,ws,rlpan)
          rlcurve0 = rlcurve0 + rlpan
!
 4600 continue
        enddo
!
        if ((ifdone .eq. 1) .and. (nchnew .eq. nch)) goto 5100
        nchnew=nch
!        print *, ijk,rlcurve,rkmax,rkmin,epsuse,2.0d0**(-51)*rkmax/rkmin*10
      enddo
!
 5100 continue
!
!       the curve should be resolved to precision eps now on
!       each interval ab(,i)
!       check the size of adjacent neighboring chunks - if off by a
!       factor of more than 2, split them as well. iterate until done.
!

      maxiter=20
      do ijk=1,maxiter
!
        nchold=nch
        ifdone=1
        do i=1,nchold
!
          i1=adjs(1,i)
          i2=adjs(2,i)

!
!       calculate chunk lengths
!
          a=ab(1,i)
          b=ab(2,i)
          rlself = 0
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,b,xs,ws,rlself)
!
          rl1=rlself
          rl2=rlself
!
          if (i1 .gt. 0) then
            a1=ab(1,i1)
            b1=ab(2,i1)
            rl1 = 0
            call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a1,b1,xs,ws,rl1)
          endif
!
          if (i2 .gt. 0) then
            a2=ab(1,i2)
            b2=ab(2,i2)
            rl2 = 0
            call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a2,b2,xs,ws,rl2)
          endif

!
!       only check if self is larger than either of adjacent blocks,
!       iterating a couple times will catch everything
!
          ifsplit=0
          sc = 2.001d0
          if (rlself .gt. sc*rl1) ifsplit=1
          if (rlself .gt. sc*rl2) ifsplit=1
          if (ifsplit .eq. 0) goto 8600
!
!       split chunk i now, and recalculate nodes, ders, etc
!
          ifdone=0
          a=ab(1,i)
          b=ab(2,i)
          ab2=(a+b)/2
!
          i1=adjs(1,i)
          i2=adjs(2,i)
!        
          adjs(1,i) = i1
          adjs(2,i) = nch+1
!
!       . . . first update nch+1
!
          nch=nch+1
          if (nch .gt. nchmax) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            ier = 4
            return
          endif
!
          adjs(1,nch) = i
          adjs(2,nch) = i2
!
!       . . . if there's an i2, update it
!
          if (i2 .gt. 0) then
            adjs(1,i2) = nch
          endif
!
          ab(1,i)=a
          ab(2,i)=ab2
!
          ab(1,nch)=ab2
          ab(2,nch)=b
!
 8600 continue
        enddo
!
        if (ifdone .ne. 0) goto 9100
      enddo
 9100 continue

!
!       go ahead and oversample by nover, updating
!       the adjacency information adjs along the way
!
      if (nover .le. 1) goto 6100
      nchtest = nch*2**(nover-1)
      if (nchtest .gt. nchmax) then
        call prinf('too many chunks in chunkfunc!*',done,0)
        ier = 4
        return
      endif

!
!
      do ijk=1,nover-1
!
        nchold=nch
        do i=1,nchold
!
          a=ab(1,i)
          b=ab(2,i)
!
!       find ab2 using newton such that 
!       len(a,ab2)=len(ab2,b)=half the chunk length
!
          rl = 0.0d0
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,b,xs,ws,rl)
!
          rlhalf=rl/2
          thresh=1.0d-8
          ifnewt=0
          ab0=(a+b)/2
!
          do iter=1,1000
!
            rl1 = 0.0d0
            call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,ab0,xs,ws,rl1)
            call funcurve(ab0,ndd,dpars,ndz,zpars,ndi,ipars,x,y, &
              dx,dy,ddx,ddy)
            dsdt=sqrt(dx**2+dy**2)
            ab1=ab0-(rl1-rlhalf)/dsdt
!
            err=rl1-rlhalf
            if (abs(err) .lt. thresh) ifnewt=ifnewt+1
!
            if (ifnewt .eq. 3) goto 6700
            ab0=ab1
          enddo
 6700 continue
!
          if (ifnewt .lt. 3) then
            call prin2('newton failed! interval not split.*',done,0)
            ier = 6
            return
          endif
!
          ab2=ab1
!
          i1=adjs(1,i)
          i2=adjs(2,i)
          adjs(2,i)=nch+1
          if (i2 .gt. 0) adjs(1,i2)=nch+1
!
          adjs(1,nch+1)=i
          adjs(2,nch+1)=i2
!
          ab(1,i)=a
          ab(2,i)=ab2
!
          nch=nch+1
          ab(1,nch)=ab2
          ab(2,nch)=b
!
        enddo
      enddo
 6100 continue
!
!
!       check the dyadic refinement at the
!       ends, first find the end segments
!
      do i =1,nch
        if (abs(ab(1,i)-ta) .lt. 1.0d-15) ileft=i
        if (abs(ab(2,i)-tb) .lt. 1.0d-15) iright=i
      enddo
!
      a1=ab(1,ileft)
      b1=ab(2,ileft)
      a2=ab(1,iright)
      b2=ab(2,iright)
!
      rlleft = 0.0d0
      rlright = 0.0d0
      call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
        a1,b1,xs,ws,rlleft)
      call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
        a2,b2,xs,ws,rlright)
!
!       . . . dyadically split the left segment
!
      if (irefinel.eq.1) then
        do ijk=1,52

!
          if (rlleft .le. rlmaxe) goto 7100
          a=ab(1,ileft)
          b=ab(2,ileft)
          ab2=(a+b)/2
!
          i1=adjs(1,ileft)
          i2=adjs(2,ileft)
          adjs(2,ileft)=nch+1
          if (i2 .gt. 0) adjs(1,i2)=nch+1
!
          nch=nch+1
          if (nch .gt. nchmax) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            ier = 4
            return
          endif
          adjs(1,nch)=ileft
          adjs(2,nch)=i2
!
          ab(1,ileft)=a
          ab(2,ileft)=ab2
!
!
          ab(1,nch)=ab2
          ab(2,nch)=b
!
          rlleft = 0.0d0
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            a,ab2,xs,ws,rlleft)
!
        enddo
      endif
 7100 continue
!
!       . . . dyadically split the right segment
!
      if(irefiner.eq.1) then
        do ijk=1,52
!
!
          if (rlright .le. rlmaxe) goto 7500
          a=ab(1,iright)
          b=ab(2,iright)
          ab2=(a+b)/2
!
          i1=adjs(1,iright)
          i2=adjs(2,iright)
          adjs(2,iright)=nch+1
          nch=nch+1
          if (nch .gt. nchmax) then
            call prinf('too many chunks in chunkfunc!*',done,0)
            return
          endif
!
          adjs(1,nch)=iright
          adjs(2,nch)=i2
!
          ab(1,iright)=a
          ab(2,iright)=ab2
!
!
          iright=nch
          ab(1,nch)=ab2
          ab(2,nch)=b
!
          rlright = 0.0d0
          call chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
            ab2,b,xs,ws,rlright)
!
        enddo
      endif
 7500 continue
!
 7700 continue


      ich = 1
      if(ifclosed.eq.0) then

        do i=1,nch
          if(adjs(1,i).lt.0) ich = i
        enddo
      endif
      

!
!      up to here, everything has been done in parameter space, [ta,tb]
!      . . . finally evaluate the k nodes on each chunk, along with 
!      derivatives and chunk lengths
!
      do i = 1, nch
!
        a=ab(1,ich)
        b=ab(2,ich)
        hs=(b-a)/2
        norders(i) = k
        ixys(i) = (i-1)*k + 1
        iptype(i) = 1
!
        do j = 1, k
          ipt = (i-1)*k + j
          t=a+(b-a)*(xs(j)+1)/2
          call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
            x,y,dx,dy,dx2,dy2)
          srcvals(1,ipt) = x
          srcvals(2,ipt) = y
          srcvals(3,ipt) = dx*hs
          srcvals(4,ipt) = dy*hs
          srcvals(5,ipt) = dx2*hs**2
          srcvals(6,ipt) = dy2*hs**2
          ds = sqrt(dx**2 + dy**2)
          srcvals(7,ipt) = dy/ds
          srcvals(8,ipt) = -dx/ds
        enddo
        call dgemm('n','t',6,k,k,alpha,srcvals(1,ixys(i)),8,u,k, &
            beta,srccoefs(1,ixys(i)),6)
        ich = adjs(2,ich)
      enddo
      ixys(nch+1) = nch*k+1
      npts = nch*k
      ier = 0
      do i=1,nch
        adjs(1,i) = i-1
        adjs(2,i) = i+1
      enddo
      if(ifclosed.eq.0) then
        adjs(1,1) = -1
        adjs(2,nch) = -1
      else
        adjs(1,1) = nch
        adjs(2,nch) = 1
      endif

!
      return
      end
!---------------------------------------------------------------------
!
!
!
!
!
  
      subroutine chunklength(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
        ta,tb,xs,ws,rl)
!
!  Using a k point gaussian quadrature, calculate the length
!  of a curve between a and b. The curve as before
!  is specificed throught a user-defined subroutine funcurve
!  whose calling expected calling sequence is
!
!  funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars,x,y,dxdt,dydt,dxdt2,dydt2)
!
!  Input arguments:  
!    - k: integer
!        number of legendre nodes to be used in gaussian quadrature 
!    - funcurve: function handle
!        function which returns x,y,dxdt,dydt,dxdt2,dydt2
!        for a given point in parameter space
!    - ndd: integer
!        length of double precision parameter array for funcurve
!    - dpars: double precision(ndd)
!        list of double precision parameters for funcurve
!    - ndz: integer
!        length of complex parameter array for funcurve
!    - dpars: double comlex(ndd)
!        list of complex parameters for funcurve
!    - ndi: integer
!        length of integer parameter array for funcurve
!    - ipars: integer(ndi)
!        list of integer parameters for funcurve
!    - ta,tb: double precision
!        find length of curve for ta\leq t \leq tb
!    - xs,ws: double precision(k)
!        order k legendre nodes and weights
!
!  Output arguments:
!    - rl: double precision
!       length of curve 

      implicit none
      integer, intent(in) :: k,ndd,ndz,ndi,ipars(ndi)
      real *8, intent(in) :: dpars(ndd),ta,tb,xs(k),ws(k)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(out) :: rl
      
      integer i
      real *8 t,dsdt,h,x,y,dx,dy,dx2,dy2

      external funcurve
      
      rl = 0
      h = (tb-ta)/2
      do i=1,k
        t = ta + (tb-ta)*(xs(i)+1)/2
        call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
              x,y,dx,dy,dx2,dy2)
        dsdt = sqrt(dx**2+dy**2)*h
        rl = rl + dsdt*ws(i)
      enddo


      return
      end
!
!
!
!
!
!
      
  
      subroutine chunkints(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
        ta,tb,xs,ws,rints,rl)
!
!  Using a k point gaussian quadrature, calculate the
!  following integrals on the curve between a,b
!  rints(1) = \int_{a}^{b} |dxdt(t)|^2 dt
!  rints(2) = \int_{a}^{b} |dydt(t)|^2 dt
!  rints(3) = \int_{a}^{b} (dsdt)^2 dt
!  rints(4) = \int_{a}^{b} (k(t)^2dsdt) dt
!  rl = \int_{a}^{b} dsdt dt
!
!
!  The curve as before
!  is specificed throught a user-defined subroutine funcurve
!  whose calling expected calling sequence is
!
!  funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars,x,y,dxdt,dydt,dxdt2,dydt2)
!
!  Input arguments:  
!    - k: integer
!        number of legendre nodes to be used in gaussian quadrature 
!    - funcurve: function handle
!        function which returns x,y,dxdt,dydt,dxdt2,dydt2
!        for a given point in parameter space
!    - ndd: integer
!        length of double precision parameter array for funcurve
!    - dpars: double precision(ndd)
!        list of double precision parameters for funcurve
!    - ndz: integer
!        length of complex parameter array for funcurve
!    - dpars: double comlex(ndd)
!        list of complex parameters for funcurve
!    - ndi: integer
!        length of integer parameter array for funcurve
!    - ipars: integer(ndi)
!        list of integer parameters for funcurve
!    - ta,tb: double precision
!        find length of curve for ta\leq t \leq tb
!    - xs,ws: double precision(k)
!        order k legendre nodes and weights
!
!  Output arguments:
!    - rints: double precision(4)
!        the integrals above
!    - rl: double precision
!       length of curve 

      implicit none
      integer, intent(in) :: k,ndd,ndz,ndi,ipars(ndi)
      real *8, intent(in) :: dpars(ndd),ta,tb,xs(k),ws(k)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(out) :: rl,rints(4)
      
      integer i
      real *8 t,dsdt,h,x,y,dx,dy,dx2,dy2,rtmp,ra
      complex *16 zd,zd2,ima
      data ima/(0.0d0,1.0d0)/

      external funcurve
      
      rl = 0
      h = (tb-ta)/2
      do i=1,4
        rints(i) = 0
      enddo
      ra = 0
      do i=1,k
        t = ta + (tb-ta)*(xs(i)+1)/2
        call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
              x,y,dx,dy,dx2,dy2)
            
        dsdt = sqrt(dx**2+dy**2)
        zd2 = dx2 + ima*dy2
        zd = dx + ima*dy
        rtmp = imag(zd2/zd)**2/dsdt
        rints(1) = rints(1) + dx**2*h*ws(i)
        rints(2) = rints(2) + dy**2*h*ws(i)
        rints(3) = rints(3) + dsdt**2*h*ws(i)
        rints(4) = rints(4) + rtmp*h*ws(i)
        rl = rl + dsdt*ws(i)*h
        ra = ra + ws(i)
      enddo


      return
      end
!
!
!
!
!
!
      
!
!
!
!
!
  
      subroutine chunkcurv(k,funcurve,ndd,dpars,ndz,zpars,ndi,ipars, &
        ta,tb,xs,ws,curv)
!
!  At the k point gaussian quadrature, calculate the curvature
!  of a curve between a and b. The curve as before
!  is specificed throught a user-defined subroutine funcurve
!  whose calling expected calling sequence is
!
!  funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars,x,y,dxdt,dydt,dxdt2,dydt2)
!
!  Input arguments:  
!    - k: integer
!        number of legendre nodes to be used in gaussian quadrature 
!    - funcurve: function handle
!        function which returns x,y,dxdt,dydt,dxdt2,dydt2
!        for a given point in parameter space
!    - ndd: integer
!        length of double precision parameter array for funcurve
!    - dpars: double precision(ndd)
!        list of double precision parameters for funcurve
!    - ndz: integer
!        length of complex parameter array for funcurve
!    - dpars: double comlex(ndd)
!        list of complex parameters for funcurve
!    - ndi: integer
!        length of integer parameter array for funcurve
!    - ipars: integer(ndi)
!        list of integer parameters for funcurve
!    - ta,tb: double precision
!        find length of curve for ta\leq t \leq tb
!    - xs,ws: double precision(k)
!        order k legendre nodes and weights
!
!  Output arguments:
!    - curv: double precision(k)
!       length of curve 

      implicit none
      integer, intent(in) :: k,ndd,ndz,ndi,ipars(ndi)
      real *8, intent(in) :: dpars(ndd),ta,tb,xs(k),ws(k)
      complex *16, intent(in) :: zpars(ndz)
      real *8, intent(out) :: curv(k)
      
      integer i
      real *8 t,dsdt,h,x,y,dx,dy,dx2,dy2

      external funcurve
      
      h = (tb-ta)/2
      do i=1,k
        t = ta + (tb-ta)*(xs(i)+1)/2
        call funcurve(t,ndd,dpars,ndz,zpars,ndi,ipars, &
              x,y,dx,dy,dx2,dy2)
        dsdt = sqrt(dx**2+dy**2)
        curv(i) = (dy2*dx-dx2*dy)/dsdt**3
      enddo


      return
      end
!
!
!
!
!
!
      subroutine interpmat_1d_hels(k,x1,x2,ximat)
      implicit real *8 (a-h,o-z)
      real *8 x1(k),x2(k),ximat(k,k),vmat_trans(k,k),rmat_trans(k,k)
      real *8 ximat_trans(k,k)
      
      vmat_trans = 1.0d0
      rmat_trans = 1.0d0
      do j=2,k
        do i=1,k
          vmat_trans(j,i) = vmat_trans(j-1,i)*x2(i)
          rmat_trans(j,i) = rmat_trans(j-1,i)*x1(i)
        enddo
      enddo

      info = 0

      call dgausselim_vec(k,vmat_trans,k,rmat_trans,info,ximat_trans, &
        dcond)
      
      do i=1,k
        do j=1,k
          ximat(i,j) = ximat_trans(j,i)
        enddo
      enddo


      return
      end
  
