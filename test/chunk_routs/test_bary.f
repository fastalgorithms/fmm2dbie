      subroutine test_bary1d(isuc)
      
      implicit real *8 (a-h,o-z)
      real *8, allocatable :: ts(:),ws(:),umat(:,:),vmat(:,:)
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),baryw(:)
      real *8, allocatable :: pols(:)
      real *8 srctmp(8)

      done = 1
      pi = atan(done)*2
      
      itype = 2
      k = 26
      allocate(ts(k),ws(k),umat(k,k),vmat(k,k))

      write(*,*) "================================"
      write(*,*) " "
      
      call legeexps(itype,k,ts,umat,vmat,ws)

      allocate(srcvals(8,k),srccoefs(6,k))
      allocate(baryw(k))
      allocate(pols(k))
      

      do i=1,k
        tt = (ts(i)+1)*pi
        srcvals(1,i) = cos(tt)
        srcvals(2,i) = sin(tt)
        srcvals(3,i) = -sin(tt)*pi
        srcvals(4,i) = cos(tt)*pi
        srcvals(5,i) = -cos(tt)*pi*pi
        srcvals(6,i) = -sin(tt)*pi*pi
        srcvals(7,i) = cos(tt)
        srcvals(8,i) = sin(tt)

        srccoefs(1,i) = 0
        srccoefs(2,i) = 0
        srccoefs(3,i) = 0
        srccoefs(4,i) = 0
        srccoefs(5,i) = 0
        srccoefs(6,i) = 0
      enddo

      do i=1,k
        do j=1,k
          do l=1,6
            srccoefs(l,i) = srccoefs(l,i) + umat(i,j)*srcvals(l,j)
          enddo
        enddo
      enddo

      
      ttest = ts(1) + 3.2d-2 
      call legepols(ttest,k-1,pols)
      
      rr = 0
      do l=1,8
        srctmp(l) = 0
      enddo


      do i=1,k
        do l=1,6
          srctmp(l) = srctmp(l) + srccoefs(l,i)*pols(i)
        enddo
      enddo
      dst = sqrt(srctmp(3)**2 + srctmp(4)**2)

      srctmp(7) = srctmp(4)/dst
      srctmp(8) = -srctmp(3)/dst

      dx = srctmp(1) - srcvals(1,1)
      dy = srctmp(2) - srcvals(2,1)
      
      r = sqrt(dx**2 + dy**2)

      rdlpex = 0.5d0*r

      err1 = abs(rdlpex-rdlp1)
      isuc = 0
      if(err1.le.1.0d-14) isuc = isuc+1
      print *, "error in rdlp through interpolation=",err1


      
      do i=1,k
        baryw(i) = 1
        do j=1,k
          if(j.ne.i) baryw(i) = baryw(i)/(ts(i)-ts(j))
        enddo
      enddo

      do l=1,8
        srctmp(l) = 0
      enddo

      ra = 0
      do i=1,k
        ra = ra + baryw(i)/(ttest-ts(i))
      enddo

      do i=1,k
        do l=1,6
          srctmp(l) = srctmp(l) + srcvals(l,i)*baryw(i)/(ttest-ts(i))
        enddo
      enddo

      do l=1,6
        srctmp(l) = srctmp(l)/ra
      enddo

      
      dst = sqrt(srctmp(3)**2 + srctmp(4)**2)

      srctmp(7) = srctmp(4)/dst
      srctmp(8) = -srctmp(3)/dst

      dx = srctmp(1) - srcvals(1,1)
      dy = srctmp(2) - srcvals(2,1)
      
      r = sqrt(dx**2 + dy**2)
      rdlp2 = (dx*srctmp(7) + dy*srctmp(8))/r

      rdlpex = 0.5d0*r

      err2 = abs(rdlpex-rdlp2)
      print *, "error in bary interpolation=",err2
      if(err2.le.1.0d-14) isuc = isuc+1

      write(*,*)  " "
      write(*,*) "=================================================="



      

      return
      end
