% Get Indices utility function
% function [xid,yid,lxid,lyid,rxid,ryid] = getIndices(f,sz)
%    f = Sampling factor, sz = size of image
%
%    xid,yid are vectors of length equal to the width and height of
%    the image respectively. Groups of pixels (i,j) that have the
%    same values of yid(j) and xid(i) form 'patches' at whose four
%    corners we have chromaticity measurements.
%
%     Moreover, in the chromatiities returned by getChrom,
%     R(yid,xid)..R(yid+1,xid+1) give the values at the four
%     corners.
%
%     lxid,lyid,rxid,ryid give the actual pixel locations in the
%     full image of these corners.
function [xidx,yidx,lxidx,lyidx,rxidx,ryidx] = getIndices(f,sz)

% Left/Top patches are 1-pixel bigger since we aren't splitting the
% 4x4 bayer block between 4 patches.
xidx = floor( ([2 2:sz(2)]-2) / f) + 1;
yidx = floor( ([2 2:sz(1)]-2) / f) + 1;

% Corresponding corners
lxidx = [min((xidx-1)*f+2,sz(2))];
rxidx = min(lxidx+f-1,sz(2));

lyidx = [min((yidx-1)*f+2,sz(1))];
ryidx = min(lyidx+f-1,sz(1));
