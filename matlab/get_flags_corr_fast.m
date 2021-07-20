function [flag,corr] = get_flags_corr_fast(chnkr,zk,targs)
%
%

  nch = chnkr.nch;
  nchp1 = nch+1;
  k = chnkr.k;
  npts = nch*k;

  srcvals = zeros(8,npts);
  hrep = reshape(repmat(chnkr.h,[1,k])',[npts,1])';
  srcvals(1:2,:) = reshape(chnkr.r,[2,npts]);
  srcvals(3:4,:) = reshape(chnkr.d,[2,npts]);
  srcvals(3,:) = srcvals(3,:).*hrep;
  srcvals(4,:) = srcvals(4,:).*hrep;
  srcvals(5:6,:) = reshape(chnkr.d2,[2,npts]);
  srcvals(5,:) = srcvals(5,:).*hrep.^2;
  srcvals(6,:) = srcvals(6,:).*hrep.^2;
  rnorms = normals(chnkr);
  rnorms = reshape(rnorms,[2,npts]);
  srcvals(7:8,:) = rnorms;


  srccoefs = zeros(6,npts);
  [rc,dc,d2c] = exps(chnkr);
  srccoefs(1:2,:) = reshape(rc,[2,npts]);
  srccoefs(3:4,:) = reshape(dc,[2,npts]);
  srccoefs(3,:) = srccoefs(3,:).*hrep;
  srccoefs(4,:) = srccoefs(4,:).*hrep;
  srccoefs(5:6,:) = reshape(d2c,[2,npts]);
  srccoefs(5,:) = srccoefs(5,:).*hrep.^2;
  srccoefs(6,:) = srccoefs(6,:).*hrep.^2;
  
  ixys = 1:k:(npts+1);
  norders = k*ones(nch,1);
  iptype = ones(nch,1);
  [ndt,nt] = size(targs);
  cms = zeros(2,nch);
  rads = zeros(nch,1);

  mex_id_ = 'get_centroid_rads2d(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], io double[xx], io double[x])';
[cms, rads] = curve_routs(mex_id_, nch, norders, ixys, iptype, npts, srccoefs, srcvals, cms, rads, 1, nch, nchp1, nch, 1, 6, npts, 8, npts, 2, nch, nch);
  rad_near = rads*1.1;

  nnz = 0;
  mex_id_ = 'findnear2dmem(i double[xx], i int[x], i double[x], i int[x], i double[xx], i int[x], io int[x])';
[nnz] = curve_routs(mex_id_, cms, nch, rad_near, ndt, targs, nt, nnz, 2, nch, 1, nch, 1, ndt, nt, 1, 1);

  nquad = nnz*k;
  wnear = complex(zeros(nquad,1));
  iind = zeros(nquad,1);
  jind = zeros(nquad,1);

  zkuse = complex(zk);
  mex_id_ = 'getnearquad_helm_comb_dir_2d_matlab(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i int[x], i int[x], i double[xx], i double[xx], i double[x], i int[x], i int[x], i dcomplex[x], io dcomplex[x], io int[x], io int[x])';
[wnear, iind, jind] = curve_routs(mex_id_, nch, norders, ixys, iptype, npts, srccoefs, srcvals, ndt, nt, targs, cms, rad_near, nnz, nquad, zkuse, wnear, iind, jind, 1, nch, nchp1, nch, 1, 6, npts, 8, npts, 1, 1, ndt, nt, 2, nch, nch, 1, 1, 1, nquad, nquad, nquad);

  flag = [];
  corr = sparse(iind,jind,wnear,nt,npts);


end  




