c
c     this file contains four user-callable subroutines
c     for finding the rectangles in a collection which
c     contain each point
c
c     complexity estimates are for m rectangles, n points
c
c     findinrectangleslow - for each point, find the
c                         rectangles which contain it
c                         straightforward O(mn) algo
c     
c     findinrectangle - for each point, find the
c                      rectangles which contain it
c                      more complicated algo,
c                      O(mlog(m)+nlog(m)) for reasonably
c                      distributed rectangles
c     
c     each of these has a _mem routine for finding memory
c     requirements  (a bit expensive, it runs the algorithm
c                     to get this)
c      
c

      subroutine findinrectangleslow(m,rects,n,ldpts,pts,iptladdr,
     1     lpr,iptrects)
c
c     for each point, find the rectangles which contain it
c
c     input
c
c     rects - (2,4,m) real *8 array of rectangle coordinates
c     pts - (ldpts,n) real *8 array of points to check
c                     pts(1:2,i) are the coordinates of point i
c     lpr - integer, the length of the iptrects array on input
c     
c     output
c
c     iptladdr - (n+1) integer array, the rectangles which
c                 contain point i are in the iptrects list in entries 
c                 iptladdr(i):iptladdr(i+1)-1. If 
c                 iptladdr(i+1) .eq. iptladdr(i), then 
c                 there point i wasn't in any rectangles
c     iptrects - integer array, the rectangles which are relevant to
c                each point. iptrects(iptladdr(i)) is the
c                first rectangle which contains point i (if 
c                point i is in any rectangles) and so on
c
      implicit real *8 (a-h,o-z)
      real *8 :: rects(2,4,m), pts(ldpts,n), rect0(2,4)
      integer :: iptladdr(n+1), iptrects(lpr)

      istart = 1
      do j = 1,n
         iptladdr(j)=istart
         do i = 1,m
            inout=-1
            call ptinrect(rects(1,1,i),pts(1,j),inout)
            if (inout .ge. 0) then
               iptrects(istart) = i
               istart = istart + 1
            endif
         enddo
      enddo
      
      iptladdr(n+1)=istart
      
      return
      end

      subroutine findinrectangleslow_mem(m,rects,n,ldpts,pts,
     1     nnz)
      
c
c     get memory requirements for findinrectangleslow
c
c     input
c
c     rects - (2,4,m) real *8 array of rectangle coordinates
c     pts - (ldpts,n) real *8 array of points to check
c                     pts(1:2,i) are the coordinates of point i
c     
c     output
c
c     nnz - integer, the length needed for the iptrects array
c            on input to findinrectangleslow
c

      implicit real *8 (a-h,o-z)
      real *8 :: rects(2,4,m), pts(ldpts,n), rect0(2,4)
      integer :: nnz

      istart = 1
      do j = 1,n
         do i = 1,m
            inout=-1
            call ptinrect(rects(1,1,i),pts(1,j),inout)
            if (inout .ge. 0) then
               istart = istart + 1
            endif
         enddo
      enddo

      nnz = istart-1
      
      return
      end

      subroutine findinrectangle_mem(m,rects,n,ldpts,pts,nnz,
     1     ier)
c
c     get memory requirements for findinrectangle
c
c     input
c
c     rects - (2,4,m) real *8 array of rectangle coordinates
c     pts - (ldpts,n) real *8 array of points to check
c                     pts(1:2,i) are the coordinates of point i
c     
c     output
c
c     nnz - integer, the length needed for the iptrects array
c               on input to findinrectangleslow
c     ier - integer, error flag
      
c
c
c     the algorithm is based on a bounding volume hierarchy
c     in which all elements are (not necessarily grid-aligned)
c     rectangles.
c
c     the original rectangles are sorted into a quad-tree based
c     on rectangle centers
c
      
      implicit real *8 (a-h,o-z)
      dimension rects(2,4,m), pts(ldpts,n)

      real *8, allocatable :: centers(:,:), rctrs(:,:), boxsize(:)
      integer, allocatable :: itree(:), irectbox(:,:), irectlist(:)
      integer, allocatable :: iptladdr(:), isortlist(:)
      real *8, allocatable :: rectboxes(:,:,:), rectsort(:,:,:)
      integer :: iptr(8)


      allocate(rctrs(2,m),iptladdr(n+1))
      call rectcenters(m,rects,rctrs)
            
      
c     build tree on rect centers

      
      ndummy=0
      idivflag=0
      ndiv=5
      nlmin=0
      nlmax=40
      ifunif=0
      iper=0

      call cpu_time(t0)
      
      call pts_tree_mem(rctrs,m,dummy,ndummy,idivflag,ndiv,nlmin,nlmax,
     1     ifunif,iper,nlevels,nboxes,ltree)
      
      allocate(itree(ltree),centers(2,nboxes),boxsize(0:nlevels))
      call pts_tree_build(rctrs,m,dummy,ndummy,idivflag,ndiv,
     1     nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,iptr,
     2     centers,boxsize)

      allocate(irectbox(2,nboxes),irectlist(m))
      call pts_tree_sort(m,rctrs,itree,ltree,nboxes,nlevels,
     1     iptr,centers,irectlist,irectbox)

      call cpu_time(t1)

c      call prin2('time building tree *',t1-t0,1)
      
      call cpu_time(t0)

      allocate(rectsort(2,4,m),isortlist(m))
      do i = 1,m
         call copyrect(rects(1,1,irectlist(i)),
     1        rectsort(1,1,i))
         isortlist(i)=i
      enddo
      
      allocate(rectboxes(2,4,nboxes))
      call buildbrh(m,rectsort,itree,nboxes,nlevels,iptr,
     1     isortlist,irectbox,rectboxes,imergtype)

      call cpu_time(t1)
      
c      call prin2('time building brh *',t1-t0,1)
      
      call cpu_time(t0)
      call findinbrh_mem(m,rectsort,itree,nboxes,nlevels,
     1     iptr,isortlist,irectbox,rectboxes,n,ldpts,pts,iptladdr,
     2     nnz,ier)

      call cpu_time(t1)

c      call prin2('time sorting into brh *',t1-t0,1)


      return
      end

      subroutine findinrectangle(m,rects,n,ldpts,pts,iptladdr,
     1     lpr,iptrects,ier)
c
c     for each point, find the rectangles which contain it
c
c     input
c
c     rects - (2,4,m) real *8 array of rectangle coordinates
c     pts - (ldpts,n) real *8 array of points to check
c                     pts(1:2,i) are the coordinates of point i
c     lpr - integer, the length of the iptrects array on input
c     
c     output
c
c     iptladdr - (n+1) integer array, the rectangles which
c                 contain point i are in the iptrects list in entries 
c                 iptladdr(i):iptladdr(i+1)-1. If 
c                 iptladdr(i+1) .eq. iptladdr(i), then 
c                 there point i wasn't in any rectangles
c     iptrects - integer array, the rectangles which are relevant to
c                each point. iptrects(iptladdr(i)) is the
c                first rectangle which contains point i (if 
c                point i is in any rectangles) and so on
c

      
c
c     the algorithm is based on a bounding volume hierarchy
c     in which all elements are (arbitrarily oriented)
c     rectangles.
c
c     the original rectangles are sorted into a quad-tree based
c     on rectangle centers
c      

      implicit real *8 (a-h,o-z)
      dimension rects(2,4,m), pts(ldpts,n)
      integer :: iptladdr(n+1), iptrects(lpr)

      real *8, allocatable :: centers(:,:), rctrs(:,:), boxsize(:)
      integer, allocatable :: itree(:), irectbox(:,:), irectlist(:)
      integer, allocatable :: isortlist(:)      
      real *8, allocatable :: rectboxes(:,:,:), rectsort(:,:,:)
      integer :: iptr(8)

      allocate(rctrs(2,m))
      call rectcenters(m,rects,rctrs)
            
      
c     build tree on rect centers

      
      ndummy=0
      idivflag=0
      ndiv=5
      nlmin=0
      nlmax=40
      ifunif=0
      iper=0

      call cpu_time(t0)
      
      call pts_tree_mem(rctrs,m,dummy,ndummy,idivflag,ndiv,nlmin,nlmax,
     1     ifunif,iper,nlevels,nboxes,ltree)
      
      allocate(itree(ltree),centers(2,nboxes),boxsize(0:nlevels))
      call pts_tree_build(rctrs,m,dummy,ndummy,idivflag,ndiv,
     1     nlmin,nlmax,ifunif,iper,nlevels,nboxes,ltree,itree,iptr,
     2     centers,boxsize)

      allocate(irectbox(2,nboxes),irectlist(m))
      call pts_tree_sort(m,rctrs,itree,ltree,nboxes,nlevels,
     1     iptr,centers,irectlist,irectbox)

      call cpu_time(t1)

c      call prin2('time building tree *',t1-t0,1)
      
      call cpu_time(t0)

      allocate(rectsort(2,4,m),isortlist(m))
      do i = 1,m
         call copyrect(rects(1,1,irectlist(i)),
     1        rectsort(1,1,i))
         isortlist(i)=i
      enddo

      allocate(rectboxes(2,4,nboxes))
      call buildbrh(m,rectsort,itree,nboxes,nlevels,iptr,
     1     isortlist,irectbox,rectboxes,imergtype)

      call cpu_time(t1)
      
c      call prin2('time building brh *',t1-t0,1)
      
      call cpu_time(t0)
      call findinbrh(m,rectsort,itree,nboxes,nlevels,
     1     iptr,isortlist,irectbox,rectboxes,n,ldpts,pts,iptladdr,
     2     lpr,iptrects,ier)

      do i = 1,(iptladdr(n+1)-1)
         iptrects(i) = irectlist(iptrects(i))
      enddo
      
      call cpu_time(t1)

c      call prin2('time sorting into brh *',t1-t0,1)
      

      return
      end

      subroutine findinbrh(m,rects,itree,nboxes,nlevels,
     1     iptr,irectlist,irectbox,rectbox,n,ldp,pts,iptladdr,
     2     lpr,iptrects,ier)
      implicit real *8 (a-h,o-z)
      real *8 :: pts(ldp,n), rectbox(2,4,nboxes), rects(2,4,m)
      integer :: iptladdr(n+1), iptrects(lpr), irectbox(2,nboxes)
      integer :: irectlist(m), iptr(8), itree(*)

      real *8 :: ptloc(2)
      
      istart = 1
      ier=0
      do j = 1,n
         iptladdr(j)=istart
         ptloc(1)=pts(1,j)
         ptloc(2)=pts(2,j)
         lpr0=lpr-istart+1
         call findinbrh0(m,rects,itree,itree(iptr(4)),itree(iptr(5)),
     1        nboxes,nlevels,iptr,irectlist,irectbox,rectbox,ptloc,
     2        nrect,lpr0,iptrects(istart),ier0)
         if (ier0 .ne. 0) then
            ier = 1
            return
         endif
         istart = istart + nrect
      enddo
      iptladdr(n+1)=istart
      
      return
      end
      
      subroutine findinbrh_mem(m,rects,itree,nboxes,nlevels,
     1     iptr,irectlist,irectbox,rectbox,n,ldp,pts,iptladdr,
     2     nnz,ier)
      implicit real *8 (a-h,o-z)
      real *8 :: pts(ldp,n), rectbox(2,4,nboxes), rects(2,4,m)
      integer :: iptladdr(n+1), irectbox(2,nboxes)
      integer :: irectlist(m), iptr(8), itree(*)

      real *8 :: ptloc(2)
      
      istart = 1
      ier=0
      do j = 1,n
         iptladdr(j)=istart
         ptloc(1)=pts(1,j)
         ptloc(2)=pts(2,j)
         lpr0=lpr-istart+1
         call findinbrh0_mem(m,rects,itree,itree(iptr(4)),
     1        itree(iptr(5)),
     1        nboxes,nlevels,iptr,irectlist,irectbox,rectbox,ptloc,
     2        nrect,ier0)
         if (ier0 .ne. 0) then
            ier = 1
            return
         endif
         istart = istart + nrect
      enddo
      iptladdr(n+1)=istart
      nnz = istart-1
      
      return
      end
      

      subroutine findinbrh0(m,rects,itree,nchild,ichild,
     1     nboxes,nlevels,iptr,irectlist,irectbox,rectbox,pt,
     2     nrect,lpr,iptrects,ier)
      
      implicit real *8 (a-h,o-z)
      real *8 :: pt(2), rectbox(2,4,nboxes), rects(2,4,m)
      integer :: iptrects(lpr), irectbox(2,nboxes)
      integer :: irectlist(m), iptr(8), itree(*), ichild(4,nboxes)
      integer :: nchild(nboxes)
c     local
      integer :: istack(200)
      real *8 :: rect0(2,4)

      ier=0
      is = 1
      istack(1)=1
      ntry=0
      np=4
      nrect=0
      ncheck=0
      
      do while(is .gt. 0 .and. ntry .le. nboxes)
         ntry=ntry+1
         ibox = istack(is)
         if (irectbox(2,ibox) .ge. irectbox(1,ibox)) then
            ncheck = ncheck+1
            call ptinrect(rectbox(1,1,ibox),pt,inout)
            if (inout .ge. 0) then
c     in box rectangle
               if (nchild(ibox) .gt. 0) then
c     has children, add em to stack
                  istack(is) = ichild(1,ibox)
                  is = is+1
                  istack(is) = ichild(2,ibox)
                  is = is+1
                  istack(is) = ichild(3,ibox)
                  is = is+1
                  istack(is) = ichild(4,ibox)
               else
c     leaf, check actual rects
                  do j = irectbox(1,ibox),irectbox(2,ibox)
                     irect = irectlist(j)
                     call ptinrect(rects(1,1,irect),pt,inout)
                     
                     if (inout .ge. 0) then
                        nrect = nrect+1
                        if (nrect .gt. lpr) then
                           ier = 1
                           exit
                        endif
                        iptrects(nrect)=irect
                     endif
                  enddo
                  is = is-1
               endif
            else
c     not in this rectangle
               is = is-1
            endif
         else
c     box doesn't even have rectangles
            is = is-1
         endif
      enddo
                        
      return
      end
      
      subroutine findinbrh0_mem(m,rects,itree,nchild,ichild,
     1     nboxes,nlevels,iptr,irectlist,irectbox,rectbox,pt,
     2     nrect,ier)
      
      implicit real *8 (a-h,o-z)
      real *8 :: pt(2), rectbox(2,4,nboxes), rects(2,4,m)
      integer :: irectbox(2,nboxes)
      integer :: irectlist(m), iptr(8), itree(*), ichild(4,nboxes)
      integer :: nchild(nboxes)
c     local
      integer :: istack(200)
      real *8 :: rect0(2,4)

      ier=0
      is = 1
      istack(1)=1
      ntry=0
      np=4
      nrect=0
      
      do while(is .gt. 0 .and. ntry .le. nboxes)
         ntry=ntry+1
         ibox = istack(is)
         if (irectbox(2,ibox) .ge. irectbox(1,ibox)) then
            call ptinrect(rectbox(1,1,ibox),pt,inout)
            if (inout .ge. 0) then
c     in box rectangle
               if (nchild(ibox) .gt. 0) then
c     has children, add em to stack
                  istack(is) = ichild(1,ibox)
                  is = is+1
                  istack(is) = ichild(2,ibox)
                  is = is+1
                  istack(is) = ichild(3,ibox)
                  is = is+1
                  istack(is) = ichild(4,ibox)
               else
c     leaf, check actual rects
                  do j = irectbox(1,ibox),irectbox(2,ibox)
                     irect = irectlist(j)
                     call ptinrect(rects(1,1,irect),pt,inout)
                     if (inout .ge. 0) then
                        nrect = nrect+1
                     endif
                  enddo
                  is = is-1
               endif
            else
c     not in this rectangle
               is = is-1
            endif
         else
c     box doesn't even have rectangles
            is = is-1
         endif
      enddo
                        
      return
      end
      

      subroutine buildbrh(m,rects,itree,nboxes,nlevels,iptr,
     1     irectlist,irectbox,rectbox,imergtype)
      implicit real *8 (a-h,o-z)
      dimension rects(2,4,m), rectbox(2,4,nboxes)
      parameter (nrectloc=10)
      real *8 s1(2), s2(2), vave(2), rectsloc(2,4,nrectloc)
      integer ic(4), itree(*), iptr(8), irectlist(m+1)
      integer irectbox(2,nboxes)
      
c     upward pass
c     
c     if a leaf, merge rectangles within
c     if not a leaf, merge rectangles in children

      do i=1,nboxes
         do j = 1,4
            rectbox(1,j,i)=0
            rectbox(2,j,i)=0
         enddo
      enddo
      
      do ilev = nlevels,0,-1
         do ibox=itree(2*ilev+1),itree(2*ilev+2)
            if(irectbox(2,ibox).ge.irectbox(1,ibox)) then
               if(itree(iptr(4)+ibox-1).gt.0) then
                  istart=iptr(5)+(ibox-1)*4-1
                  ic(1) = itree(istart+1)
                  ic(2) = itree(istart+2)
                  ic(3) = itree(istart+3)
                  ic(4) = itree(istart+4)

                  nrect=0
                  do j = 1,4
                     if (irectbox(2,ic(j)).ge.irectbox(1,ic(j))) then
                        nrect=nrect+1
                        call copyrect(rectbox(1,1,ic(j)),
     1                       rectsloc(1,1,nrect))
                     endif
                  enddo
                  call boundingrect(nrect,rectsloc,rectbox(1,1,ibox))
               else
                  nrect=0
                  do j = irectbox(1,ibox),irectbox(2,ibox)
                     nrect=nrect+1
                     irect = irectlist(j)
                     call copyrect(rects(1,1,irect),
     1                    rectsloc(1,1,nrect))
                  enddo
                  call boundingrect(nrect,rectsloc,rectbox(1,1,ibox))
                     
               endif
            endif
         enddo
      enddo

      return
      end

      subroutine boundingrect(m,rects,brect)
      implicit real *8 (a-h,o-z)
      dimension rects(2,4,m), brect(2,4)
      real *8 :: recttemp(2,4)

      if (m .lt. 1) return

      call copyrect(rects,brect)
      
      do i = 2,m
         call boundingrect2(brect,rects(1,1,i),recttemp)
         call copyrect(recttemp,brect)
      enddo
      
      return
      end

      subroutine copyrect(rectin,rectout)
      implicit real *8 (a-h,o-z)
      real *8 :: rectin(2,4), rectout(2,4)
      
      do j = 1,4
         do i = 1,2
            rectout(i,j)=rectin(i,j)
         enddo
      enddo
      return
      end
      
      subroutine boundingrect2(rect1,rect2,brect)
c
c     This routine computes the optimal bounding rectangle
c     for two rectangles. The method is brutish but
c     the code is pretty lean.
c
c
c     The method is based on the "rotating calipers" idea:
c     the optimal rectangle will have to have a side
c     which is parallel to an edge of the convex
c     hull of the two rectangles.
c      
      implicit real *8 (a-h,o-z)
      dimension rect1(2,4), rect2(2,4), brect(2,4)
      real *8 :: vs(2,20), vperp(2,20), coords1(8)
      real *8 :: coords2(8)

      ii = 0
      do i = 1,4
         x = rect2(1,i)
         y = rect2(2,i)
         do j = 1,4
            ii = ii+1
            vs(1,ii) = rect1(1,j)-x
            vs(2,ii) = rect1(2,j)-y
         enddo
      enddo

      vs(1,17) = rect1(1,1)-rect1(1,2)
      vs(2,17) = rect1(2,1)-rect1(2,2)
      vs(1,18) = rect1(1,3)-rect1(1,2)
      vs(2,18) = rect1(2,3)-rect1(2,2)
      vs(1,19) = rect2(1,1)-rect2(1,2)
      vs(2,19) = rect2(2,1)-rect2(2,2)
      vs(1,20) = rect2(1,3)-rect2(1,2)
      vs(2,20) = rect2(2,3)-rect2(2,2)

      amin=huge(1d0)
      do j = 1,20
         v1 = vs(1,j)
         v2 = vs(2,j)
         vp1 = -v2
         vp2 = v1
         do i = 1,4
            coords1(i)=v1*rect1(1,i)+v2*rect1(2,i)
            coords2(i)=vp1*rect1(1,i)+vp2*rect1(2,i)
            coords1(i+4)=v1*rect2(1,i)+v2*rect2(2,i)
            coords2(i+4)=vp1*rect2(1,i)+vp2*rect2(2,i)
         enddo
         dbig1=-huge(1d0)
         dsml1=huge(1d0)
         dbig2=-huge(1d0)
         dsml2=huge(1d0)
         do i = 1,8
            dbig1=max(dbig1,coords1(i))
            dsml1=min(dsml1,coords1(i))
            dbig2=max(dbig2,coords2(i))
            dsml2=min(dsml2,coords2(i))
         enddo
         vv = v1*v1+v2*v2
         aj = abs( (dbig1-dsml1)*(dbig2-dsml2)/vv )
         if (aj .lt. amin) then
            amin = aj
            dbig1=dbig1/vv
            dsml1=dsml1/vv
            dbig2=dbig2/vv
            dsml2=dsml2/vv
            brect(1,1) = v1*dbig1+vp1*dbig2
            brect(2,1) = v2*dbig1+vp2*dbig2
            brect(1,2) = v1*dbig1+vp1*dsml2
            brect(2,2) = v2*dbig1+vp2*dsml2
            brect(1,3) = v1*dsml1+vp1*dsml2
            brect(2,3) = v2*dsml1+vp2*dsml2
            brect(1,4) = v1*dsml1+vp1*dbig2
            brect(2,4) = v2*dsml1+vp2*dbig2
         endif
      enddo
      
      return
      end


      subroutine rectcenters(m,rects,ctrs)
      implicit real *8 (a-h,o-z)
      dimension rects(2,4,m), ctrs(2,m)


      do i = 1,m
         ctrs(1,i)=0
         ctrs(2,i)=0
         do j = 1,4
            ctrs(1,i)=ctrs(1,i)+rects(1,j,i)
            ctrs(2,i)=ctrs(2,i)+rects(2,j,i)
         enddo
         ctrs(1,i)=ctrs(1,i)/4
         ctrs(2,i)=ctrs(2,i)/4
      enddo

      return
      end

      subroutine rect2file(m,rects,iw)
      implicit real *8 (a-h,o-z)
      real *8 :: rects(2,4,m)

 1000 format(e41.35)

      do j = 1,m
         do i = 1,4
            write(iw,1000) rects(1,i,j)
            write(iw,1000) rects(2,i,j)
         enddo
      enddo

      return
      end
      

      subroutine ptinrect(rect,pt,inout)
      implicit real *8 (a-h,o-z)
      real *8 :: rect(2,4), pt(2)

      x = rect(1,2)
      y = rect(2,2)
      vx = rect(1,1)-x
      vy = rect(2,1)-y
      wx = rect(1,3)-x
      wy = rect(2,3)-y

      ux = pt(1)-x
      uy = pt(2)-y

      d1=ux*vx+uy*vy
      d2=ux*wx+uy*wy

      vv=vx*vx+vy*vy
      ww=wx*wx+wy*wy

      inout = -1
      if ((d1 .ge. 0) .and. (d1 .le. vv) .and. (d2 .ge. 0) .and.
     1     (d2 .le. ww)) inout=1

      return
      end
      
