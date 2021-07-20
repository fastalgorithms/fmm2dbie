function in = chunkerinterior_fmm(chnkr,pts)
%
%CHUNKERINTERIOR_FMM returns an array indicating whether each point in
% pts is inside the domain. Assumes the domain is closed.
%
% Syntax: in = chunkerinterior(chnkr,pts)
%
% Input:
%   chnkr - chunker object describing geometry
%   pts - (chnkr.dim,:) array of points to test
%
% Output:
%   in - logical array, if in(i) is true, then pts(:,i) is inside the
%       domain
%
% Examples:
%   chnkr = chunkerfunc(@(t) starfish(t));
%   pts = 2*randn(2,100);
%   in = chunkerinterior_fmm(chnkr,pts);
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
  [ndt,nt] = size(pts);
  inflag = zeros(nt,1);
  mex_id_ = 'chunk_interior(i int[x], i int[x], i int[x], i int[x], i int[x], i double[xx], i double[xx], i double[xx], i int[x], i int[x], io int[x])';
[inflag] = curve_routs(mex_id_, nch, norders, ixys, iptype, npts, srcvals, srccoefs, pts, ndt, nt, inflag, 1, nch, nchp1, nch, 1, 8, npts, 6, npts, ndt, nt, 1, 1, nt);
  in = abs(inflag-1) < abs(inflag+1);

end  
