      implicit real *8 (a-h,o-z)
      real *8, allocatable :: rects(:,:,:), pts(:,:)
      integer, allocatable :: iptladdr(:), iptrects(:)
      integer, allocatable :: iptladdr2(:), iptrects2(:)

      call prini(6,13)
      
      iseed=1234
      dburn = hkrand(iseed)
      pi = 4*atan(1d0)

      ifcomparetoslow=0
      m = 10000
      n = 1000000
      allocate(rects(2,4,m),pts(2,n))

      do j = 1,n
         x = -1+2*hkrand(0)
         y = -1+2*hkrand(0)
         r = sqrt(x**2+y**2)
         shift = 0.9d0+0.2d0*hkrand(0)
         pts(1,j)=x*shift/r
         pts(2,j)=y*shift/r
         pts(1,j)=x
         pts(2,j)=y
      enddo

      do i = 1,m
         th = i*2*pi/m
         x = cos(th)
         y = sin(th)
         tx = -y
         ty = x

         w = (1+hkrand())*2*pi/m
         h = (1+hkrand())*0.2d0

         rects(1,1,i) = x+w*tx+h*x
         rects(2,1,i) = y+w*ty+h*y
         rects(1,2,i) = x+w*tx-h*x
         rects(2,2,i) = y+w*ty-h*y
         rects(1,3,i) = x-w*tx-h*x
         rects(2,3,i) = y-w*ty-h*y
         rects(1,4,i) = x-w*tx+h*x
         rects(2,4,i) = y-w*ty+h*y
      enddo

      l2=2

      if (ifcomparetoslow .eq. 1) then 
         call findinrectangleslow_mem(m,rects,n,l2,pts,nnz)

         call prinf('nnz *',nnz,1)

         allocate(iptladdr2(n+1),iptrects2(nnz))

         call cpu_time(t0)
         call findinrectangleslow(m,rects,n,l2,pts,iptladdr2,
     1        nnz,iptrects2)
         call cpu_time(t1)
         call prin2('time, slow rout *',t1-t0,1)
      endif

      allocate(iptladdr(n+1))
      call cpu_time(t0)
      call findinrectangle_mem(m,rects,n,l2,pts,nnz,
     1     ier)
      call cpu_time(t1)

      call prin2('time, fast rout mem*',t1-t0,1)

      write(*,*) nnz 
      
      call cpu_time(t0)
      allocate(iptrects(nnz))
      call findinrectangle(m,rects,n,l2,pts,iptladdr,nnz,iptrects,ier)
      call cpu_time(t1)

      call prin2('time, fast rout *',t1-t0,1)
      
      write(*,*) iptladdr(n+1)-1

      if (ifcomparetoslow .eq. 1) then
         call test_compareladdrs(n,iptladdr,iptladdr2,iptrects,
     1        iptrects2,ier)
         
         call prinf('ier compare *',ier,1)
      endif
      
      stop
      end


      subroutine test_compareladdrs(n,iptladdr,iptladdr2,iptrects,
     1     iptrects2,ier)
      implicit real *8 (a-h,o-z)
      integer :: iptladdr(n+1),iptladdr2(n+1),iptrects(*),iptrects2(*)

      ier = 0
      
      do i = 1,n+1
         if (iptladdr(i).ne. iptladdr2(i)) ier=1
      enddo

      ntest=0
      do i = 1,n
         do j = iptladdr(i),(iptladdr(i+1)-1)
            ntest=ntest+1
            irect = iptrects(j)
            ifound = 0
            do k = iptladdr2(i),(iptladdr2(i+1)-1)
               irect2 = iptrects2(k)
               if (irect .eq. irect2) ifound=1
            enddo
            if (ifound .eq. 0) ier = 2
         enddo
      enddo

      return
      end
      
