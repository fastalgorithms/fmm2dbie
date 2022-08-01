      implicit real *8 (a-h,o-z)
      real *8, allocatable :: srcvals(:,:),rdlpcoefs(:),pols(:)
      real *8, allocatable :: rspcoefs(:)
      real *8, allocatable :: ts(:),wts(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: xnew(:),ynew(:),dxdtnew(:),dydtnew(:)
      real *8, allocatable :: srcvals_new(:,:)
      real *8, allocatable :: xmat(:,:),rhs(:,:),soln(:,:)
      real *8, allocatable :: xint(:,:),d2xc(:),d2yc(:),dxc(:),dyc(:)
      real *8, allocatable :: xintnew(:,:,:)
      real *8, allocatable :: rdlpcoefs_all(:,:),rspcoefs_all(:,:)
      
      call prini(6,13)
      done = 1.0d0
      pi = atan(done)*4
      

      k = 16
      allocate(ts(k),wts(k),umat(k,k),vmat(k,k))
      allocate(xnew(k),ynew(k),dxdtnew(k),dydtnew(k))
      itype = 2
      call legeexps(itype,k,ts,umat,vmat,wts)

      h = 0.0000001d0
      tstart = pi/8
      allocate(srcvals(8,k),rdlpcoefs(k),pols(k),rspcoefs(k))
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
      print *, "err1=",abs(dd+0.5d0)*h
      print *, "err2=",abs(rr+0.5d0)*h

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
      rr = 0
      do i=1,k
        rr = rr + rdlpcoefs(i)*pols(i)
      enddo
      print *, "err3=",abs(rr+0.5d0)*h

      allocate(xmat(k,k-2),rhs(k,2),soln(k-2,2))
      do i=1,k
        call legepols(ts(i),k-3,pols)
        do j=1,k-2
          xmat(i,j) = pols(j)*(ts(i)-ts(ipt))**2
        enddo
      enddo

      do i=1,k
        rhs(i,1)=(srcvals(1,ipt)-srcvals(1,i))*srcvals(7,i)+
     1      (srcvals(2,ipt)-srcvals(2,i))*srcvals(8,i)
        rhs(i,2)=(srcvals(1,i)-srcvals(1,ipt))**2 + 
     1      (srcvals(2,i)-srcvals(2,ipt))**2
      enddo


      eps = 1.0d-14
      call dleastsq(k,k-2,xmat,2,rhs,eps,info,soln,irank)
      call prinf('irank=*',irank,1)
      call prin2('soln=*',soln,2*(k-2))
      
      
      call legepols(tloc,k-1,pols)

      rnum = 0
      rden = 0
      do i=1,k-2
        rnum = rnum + soln(i,1)*pols(i)
        rden = rden + soln(i,2)*pols(i)
      enddo
      rfac = rnum/rden
      call prin2_long('rfac=*',rfac,1)
      print *, "err4=",abs(rfac+0.5d0)*h

      do i=1,k
        rhs(i,1)=(srcvals_new(1,ipt)-srcvals_new(1,i))*srcvals_new(7,i)+
     1      (srcvals_new(2,ipt)-srcvals_new(2,i))*srcvals_new(8,i)
        rhs(i,2)=(srcvals_new(1,i)-srcvals_new(1,ipt))**2 + 
     1      (srcvals_new(2,i)-srcvals_new(2,ipt))**2
      enddo


      eps = 1.0d-14
      call dleastsq(k,k-2,xmat,2,rhs,eps,info,soln,irank)
      call prinf('irank=*',irank,1)
      call prin2('soln=*',soln,2*(k-2))
      
      
      call legepols(tloc,k-1,pols)

      rnum = 0
      rden = 0
      do i=1,k-2
        rnum = rnum + soln(i,1)*pols(i)
        rden = rden + soln(i,2)*pols(i)
      enddo
      rfac = rnum/rden
      call prin2_long('rfac=*',rfac,1)
      print *, "err5=",abs(rfac+0.5d0)*h

      allocate(xint(k,k),xintnew(k,k,k))
      call legeinmt_allnodes(k,xintnew)
      xint = 0

      do inode=1,k
         tstart = ts(ipt)
         tend = ts(inode)
         hh = tend - tstart
         do l=1,k
           tuse = (tend+tstart)/2 + hh/2*ts(l) 
           call legepols(tuse,k-1,pols)
           xint(1:k,inode) = xint(1:k,inode) + pols(1:k)*wts(l)*hh/2
         enddo
      enddo

c
c  
c
      do i=1,k
        srcvals_new(5,i) = srcvals(5,i)*srcvals(7,ipt) + 
     1      srcvals(6,i)*srcvals(8,ipt)
        srcvals_new(6,i) = -srcvals(5,i)*srcvals(8,ipt) + 
     1      srcvals(6,i)*srcvals(7,ipt)
      enddo
      allocate(d2xc(k),d2yc(k))
      alpha = 1.0d0
      beta = 0.0d0
      
      call dgemv('n',k,k,alpha,umat,k,srcvals_new(5,1:k),1,beta,
     1   d2xc,1)

      call dgemv('n',k,k,alpha,umat,k,srcvals_new(6,1:k),1,beta,
     1   d2yc,1)
      
      srcvals_new(1:4,1:k) = 0
      do inode=1,k
        do l=1,k
          srcvals_new(3,inode) = srcvals_new(3,inode) + 
     1       xintnew(l,inode,ipt)*d2xc(l)
          srcvals_new(4,inode) = srcvals_new(4,inode) + 
     1       xintnew(l,inode,ipt)*d2yc(l)
        enddo
      enddo

      allocate(dxc(k),dyc(k))
      call dgemv('n',k,k,alpha,umat,k,srcvals_new(3,1:k),1,beta,
     1   dxc,1)

      call dgemv('n',k,k,alpha,umat,k,srcvals_new(4,1:k),1,beta,
     1   dyc,1)
      
      do inode=1,k
        do l=1,k
          srcvals_new(1,inode) = srcvals_new(1,inode) + 
     1       xintnew(l,inode,ipt)*dxc(l)
          srcvals_new(2,inode) = srcvals_new(2,inode) + 
     1       xintnew(l,inode,ipt)*dyc(l)  
        enddo
      enddo

      do inode=1,k
        srcvals_new(4,inode) = srcvals_new(4,inode)+
     1     sqrt(srcvals(3,ipt)**2 + srcvals(4,ipt)**2)
        srcvals_new(2,inode) = srcvals_new(2,inode) + 
     1     sqrt(srcvals(3,ipt)**2 + srcvals(4,ipt)**2)*
     2      (ts(inode)-ts(ipt))
        ds = sqrt(srcvals_new(3,inode)**2 + srcvals_new(4,inode)**2)
        srcvals_new(7,inode) = srcvals_new(4,inode)/ds
        srcvals_new(8,inode) = - srcvals_new(3,inode)/ds
      enddo
      call chunk_to_lapdlpcoef_lege(k,srcvals_new,ipt,umat,rdlpcoefs)
      call prin2_long('rdlpcoefs=*',rdlpcoefs,k)
      rr = 0
      do i=1,k
        rr = rr + rdlpcoefs(i)*pols(i)
      enddo
      print *, "h=",h
      print *, "err7=",abs(rr+0.5d0)*h

      call chunk_to_lapspcoef_lege(k,srcvals_new,ipt,umat,rspcoefs)
      call prin2_long('rspcoefs=*',rspcoefs,k)
      rr = 0
      do i=1,k
        rr = rr + rspcoefs(i)*pols(i)
      enddo
      print *, "h=",h
      print *, "err8=",abs(rr-0.5d0)*h


      allocate(rdlpcoefs_all(k,k),rspcoefs_all(k,k))
      call chunk_to_ldlp_sp_xint(k,ts,srcvals,umat,xintnew,
     1   rdlpcoefs_all,rspcoefs_all)
      rr = 0
      do i=1,k
        rr = rr + rdlpcoefs_all(i,ipt)*pols(i)
      enddo
      print *, "h=",h
      print *, "err9=",abs(rr+0.5d0)

      rr = 0
      do i=1,k
        rr = rr + rspcoefs_all(i,ipt)*pols(i)
      enddo
      print *, "h=",h
      print *, "err10=",abs(rr-0.5d0)


      stop
      end


