      implicit real *8 (a-h,o-z) 
      real *8, allocatable :: srcvals(:,:),srccoefs(:,:),targs(:,:)
      real *8, allocatable :: wts(:)
      real *8, allocatable :: cms(:,:),rads(:),rad_near(:)
      integer, allocatable :: adjs(:,:)
      real *8 errs(6),ts(2)
      real *8, allocatable :: rfacs(:,:)
      character *100 fname
      integer ipars(2)
      integer, allocatable :: row_ptr(:),col_ind(:)
      integer, allocatable :: iquad(:)
      real *8, allocatable :: srcover(:,:),wover(:)
      real *8, allocatable :: uval(:),dudnval(:)
      real *8, allocatable :: slp_near(:),dlp_near(:)
      real *8, allocatable :: pot(:),potslp(:),potdlp(:)
      real *8, allocatable :: potslp2(:)

      integer, allocatable :: norders(:),ixys(:),iptype(:)
      integer, allocatable :: ixyso(:),novers(:)

      integer, allocatable :: ich_id(:),inode_id(:)
      real *8, allocatable :: ts_targ(:)
      real *8, allocatable :: ab(:,:), rects(:,:,:)
      real *8 xyz_out(2),xyz_in(2)
      real *8, allocatable :: sigma(:)
      real *8 dparsgeo(2), dpars(2)
      complex *16 zpars
      
      external fstarn_simple



      call prini(6,13)

      done = 1
      pi = atan(done)*4

      nchmax = 20000
      k = 16
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
      
      dparsgeo(1) = 0.1d0
      dparsgeo(2) = 0.05d0

      ipars(1) = 3
      nover = 1
      nch = 0
      ier = 0
      eps = 1.0d-6
      call chunkfunc_guru(eps,rlmax,ifclosed,irefinel,irefiner,rlmaxe,
     1  ta,tb,fstarn_simple,ndd,dparsgeo,ndz,zpars,ndi,ipars,nover,
     2  k,nchmax,nch,norders,ixys,iptype,npts,srcvals,srccoefs,ab,adjs,
     3  ier)

      call prinf('ier *',ier,1)
      call prinf('nch *',nch,1)
      call prinf('npts *',npts,1)

      dpars(1) = 1
      dpars(2) = 0
      
      xyz_out(1) = 3.17d0
      xyz_out(2) = -0.03d0

      xyz_in(1) = 0.17d0
      xyz_in(2) = 0.23d0

      allocate(targs(2,npts))
      allocate(wts(npts))
      call get_qwts2d(nch,norders,ixys,iptype,npts,srcvals,wts)


      allocate(cms(2,nch),rads(nch),rad_near(nch))
      allocate(pot(npts),potslp(npts),potdlp(npts))

      call get_centroid_rads2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,cms,rads)

      allocate(sigma(npts),uval(npts),dudnval(npts))

      do i=1,npts
         call l2d_slp(xyz_out,2,srcvals(1,i),1,dpars,0,zpars,0,
     1        ipars,uval(i))
         call l2d_sprime(xyz_out,8,srcvals(1,i),1,dpars,0,zpars,0,
     1        ipars,dudnval(i))
      enddo

      ndtarg = 2
     
      do i=1,npts
        targs(1,i) = srcvals(1,i)
        targs(2,i) = srcvals(2,i)
      enddo

      allocate(ich_id(npts),ts_targ(npts))
      do i=1,npts
        ich_id(i) = -1
        ts_targ(i) = 0
      enddo

      call get_chunk_id_ts(nch,norders,ixys,iptype,npts, 
     1         ich_id,ts_targ)
 
c
c    find near field
c

      rho = 2d0
      npolyfac=2

      allocate(novers(nch),ixyso(nch+1))
      
      call ellipse_nearfield2d_getnovers(eps,rho,npolyfac,
     1     nch,norders,ising,novers,ier)

      ixyso(1)=1
      do i = 1,nch
         ixyso(i+1)=ixyso(i)+novers(i)
      enddo
      npts_over=ixyso(nch+1)-1

      allocate(rects(2,4,nch))
      call ellipse_nearfield2d_definerects(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,rho,rects)

      call findinrectangle_mem(nch,rects,npts,ndtarg,targs,nnz,
     1     ier)

      allocate(row_ptr(npts+1),col_ind(nnz))

      call findinrectangle(nch,rects,npts,ndtarg,targs,row_ptr,
     1     nnz,col_ind,ier)

      allocate(iquad(nnz+1)) 
      call get_iquad_rsc2d(nch,ixys,npts,nnz,row_ptr,col_ind,
     1     iquad)

      nquad = iquad(nnz+1)-1
      allocate(slp_near(nquad),dlp_near(nquad))


      ndtarg = 2

      eps = 0.50001d-8

      ikerorder = -1


      npts_over = ixyso(nch+1)-1

      allocate(srcover(8,npts_over),wover(npts_over))

          
      call oversample_geom2d(nch,norders,ixys,iptype,npts, 
     1     srccoefs,srcvals,novers,ixyso,npts_over,srcover)
      
      call get_qwts2d(nch,novers,ixyso,iptype,npts_over,
     1     srcover,wover)


      do i=1,nquad
        slp_near(i) = 0
        dlp_near(i) = 0
      enddo



      call cpu_time(t1)

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0

      iquadtype = 1

      call getoncurvequad_lap_comb_dir_2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,
     2     eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     3     iquad,nquad,slp_near)
      
c      call getnearquad_lap_comb_dir_2d(nch,norders,
c     1     ixys,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
c     1     ich_id,ts_targ,eps,dpars,iquadtype,nnz,row_ptr,col_ind,
c     1     iquad,nquad,slp_near)

      dpars(1) = 0.0d0
      dpars(2) = 1.0d0

      call getoncurvequad_lap_comb_dir_2d(nch,norders,
     1     ixys,iptype,npts,srccoefs,srcvals,adjs,
     2     eps,dpars,iquadtype,nnz,row_ptr,col_ind,
     3     iquad,nquad,dlp_near)
      
c      call getnearquad_lap_comb_dir_2d(nch,norders,
c     1      ixys,iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
c     1      ich_id,ts_targ,eps,dpars,iquadtype,
c     1      nnz,row_ptr,col_ind,iquad,nquad,dlp_near)
      
      call cpu_time(t2)
      tquadgen = t2-t1
      
      ifinout = 0     

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0


      call cpu_time(t1)

      call lpcomp_lap_comb_dir_addsub_2d(nch,norders,ixys,
     1     iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2     eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,slp_near,
     3     dudnval,novers,npts_over,ixyso,srcover,wover,potslp)


      dpars(1) = 0.0d0
      dpars(2) = 1.0d0


      call lpcomp_lap_comb_dir_addsub_2d(nch,norders,ixys,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,
     2  eps,dpars,nnz,row_ptr,col_ind,iquad,nquad,dlp_near,
     3  uval,novers,npts_over,ixyso,srcover,wover,potdlp)


      call cpu_time(t2)
      tlpcomp = t2-t1


c
c
c      compute error
c
      errl2 = 0
      rl2 = 0
      do i=1,npts
        pot(i) = (potslp(i) - potdlp(i))*2
        errl2 = errl2 + abs(uval(i)-pot(i))**2*wts(i)
        rl2 = rl2 + abs(uval(i))**2*wts(i)

        if(i.lt.5) then
c           write(*,'(3(2x,e11.5))') potslp(i),
c     1          potdlp(i),uval(i)
        endif
        write(33,'(3(2x,e11.5))') potslp(i),
     1     potdlp(i),uval(i)
      enddo


      err = sqrt(errl2/rl2)

      call prin2('error in greens identity=*',err,1)

      i1 = 0
      if(err.lt.1.0d-2) i1 = 1

      allocate(potslp2(npts))

      dpars(1) = 1.0d0
      dpars(2) = 0.0d0
      
      call lpcomp_lap_comb_dir_2d(nch,norders,ixys,
     1  iptype,npts,srccoefs,srcvals,ndtarg,npts,targs,ich_id,
     2  ts_targ,eps,dpars,dudnval,potslp2)


      errl2 = 0
      rl2 = 0
      do i=1,npts

        errl2 = errl2 + abs(potslp(i)-potslp2(i))**2*wts(i)
        rl2 = rl2 + abs(potslp(i))**2*wts(i) 
      enddo
      errl2 = sqrt(errl2/rl2)

      call prin2('error in simpler calling interface for lp eval=*',
     1   errl2,1)

      i2 = 0
      if(errl2.lt.1.0d-12) i2 = 1
      ntests = 2

      nsuccess = i1+i2

      open(unit=33,file='../../print_testres.txt',access='append')
      write(33,'(a,i1,a,i1,a)') 'Successfully completed ',nsuccess,
     1  ' out of ',ntests,' in lap_wrappers testing suite'
      close(33)
c
      stop
      end



c
c
c
c
c
      subroutine fstarn_simple(t,ndd,dpars,ndz,zpars,ndi,ipars,x,
     1  y,dxdt,dydt,dxdt2,dydt2)
c
c   curve parametrization for simple star shaped domain
c
c       r(t) = r1 + r2*cos(n*t)
c       x(t) = r(t) cos(t)
c       y(t) = r(t) sin(t)
c   
c   where r1 = dpars(1), r2 = dpars(2) and n = ipars(1)
c     
c


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
      zd2 = ima*zd - (n+0.0d0)**2*r2*cos((n+0.0d0)*t)*exp(ima*t) -
     1    ima*(n+0.0d0)*r2*sin((n+0.0d0)*t)*exp(ima*t)
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
      end



