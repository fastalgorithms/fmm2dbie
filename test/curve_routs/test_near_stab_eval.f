      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcvals(:,:),rdlpcoefs(:),pols(:)
      real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: xnew(:),ynew(:),dxdtnew(:),dydtnew(:)
      real *8, allocatable :: srcvals_new(:,:)
      
      call prini(6,13)
      done = 1.0d0
      pi = atan(done)*4
      

      k = 16
      allocate(ts(k),wts(k),umat(k,k),vmat(k,k))
      allocate(xnew(k),ynew(k),dxdtnew(k),dydtnew(k))
      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)

      h = 0.001d0
      tstart = pi/8
      allocate(srcvals(8,k),rdlpcoefs(k),pols(k))
      allocate(srcvals_new(8,k))
      tend = tstart + h
      do i=1,k
        tuse = tstart + (ts(i)+1)/2*h
        srcvals(1,i) = cos(tuse)
        srcvals(2,i) = sin(tuse)
        srcvals(3,i) = -sin(tuse)*h/2
        srcvals(4,i) = cos(tuse)*h/2
        srcvals(5,i) = -cos(tuse)*h*h/2/2
        srcvals(6,i) = -sin(tuse)*h*h/2/2
        srcvals(7,i) = cos(tuse)
        srcvals(8,i) = sin(tuse)
      enddo

      ipt = 3
      t1 = tstart + (ts(ipt-1)+1)/2*h
      t2 = tstart + (ts(ipt+1)+1)/2*h

      tpt = t1 + hkrand(0)*(t2-t1)
      tloc = (tpt-tstart)/h*2 -1.0d0
      x = srcvals(1,ipt)-cos(tpt)
      y = srcvals(2,ipt)-sin(tpt)
      rnx = cos(tpt)
      rny = sin(tpt)
      r2 = x**2 + y**2

      
      dd = (x*rnx + y*rny)/r2
      print *, "dd=",dd
      print *, "ipt=",ipt

      call chunk_to_lapdlpcoef_lege(k,srcvals,ipt,umat,rdlpcoefs)
      call prin2('rdlpcoefs=*',rdlpcoefs,k)
      print *, "tloc=",tloc
      call legepols(tloc,k-1,pols)
      rr = 0
      do i=1,k
        rr = rr + rdlpcoefs(i)*pols(i)
      enddo
      call prin2_long('rr=*',rr,1)
      print *, "err1=",abs(dd+0.5d0)
      print *, "err2=",abs(rr+0.5d0)

c
c
c
      do i=1,k
        xnew(i) = srcvals(1,i)*srcvals(7,ipt) +
     1       srcvals(2,i)*srcvals(8,ipt)     
        ynew(i) = -srcvals(1,i)*srcvals(8,ipt) + 
     1      srcvals(2,i)*srcvals(7,ipt)
        dxdtnew(i) = srcvals(3,i)*srcvals(7,ipt) + 
     1     srcvals(4,i)*srcvals(8,ipt)
        dydtnew(i) = -srcvals(3,i)*srcvals(8,ipt) + 
     1     srcvals(4,i)*srcvals(7,ipt)
        if(i.eq.ipt) dxdtnew(i) = 0
        if(i.eq.ipt) dydtnew(i) = sqrt(srcvals(3,ipt)**2 + 
     1      srcvals(4,ipt)**2)
        srcvals_new(1,i) = xnew(i)
        srcvals_new(2,i) = ynew(i)
        srcvals_new(3,i) = dxdtnew(i)
        srcvals_new(4,i) = dydtnew(i)
        srcvals_new(5,i) = srcvals(5,i)*srcvals(7,ipt) + 
     1      srcvals(6,i)*srcvals(8,ipt)
        srcvals_new(6,i) = -srcvals(5,i)*srcvals(8,ipt) + 
     1      srcvals(6,i)*srcvals(7,ipt)
        ds = sqrt(dxdtnew(i)**2 + dydtnew(i)**2)
        srcvals_new(7,i) = dydtnew(i)/ds
        srcvals_new(8,i) = -dxdtnew(i)/ds
      enddo
      call chunk_to_lapdlpcoef_lege(k,srcvals_new,ipt,umat,rdlpcoefs)
      call prin2_long('rdlpcoefs=*',rdlpcoefs,k)
      

      


      stop
      end
