% Material Assignment / "Soft" Super-pixels  
% Assign affinities of pixels in a patch to the four corners with
% chromaticity measurements.
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function [alpha,albl] = Msegment(emap,xidx,yidx,lxidx,lyidx,rxidx,ryidx)

MAXITER=200; TOL=1e-8;

almask = ones(size(emap{1}));
almask(lyidx,lxidx) = 0; almask(lyidx,rxidx) = 0;
almask(ryidx,lxidx) = 0; almask(ryidx,rxidx) = 0;

alpha0 = zeros([size(emap{1}) 4]);
alpha0(lyidx,lxidx,1) = 1;alpha0(lyidx,rxidx,2) = 1;
alpha0(ryidx,lxidx,3) = 1;alpha0(ryidx,rxidx,4) = 1;

% Permute for label consistency
evenx = find(mod(xidx,2) == 0);
eveny = find(mod(yidx,2) == 0);
alpha0(:,evenx,:) = alpha0(:,evenx,[2 1 4 3]);
alpha0(eveny,:,:) = alpha0(eveny,:,[3 4 1 2]);
%%%

alpha = 1/4*repmat(almask,[1 1 4]) + alpha0;
alpha(:,:,1:3) = conjSolve(alpha(:,:,1:3),emap,almask,MAXITER,TOL);
alpha(:,:,4) = max(0,1-sum(alpha(:,:,1:3),3));

if nargout == 2
  [tmp,albl] = max(alpha,[],3); clear tmp
end;


% Permute back
alpha(:,evenx,:) = alpha(:,evenx,[2 1 4 3]);
alpha(eveny,:,:) = alpha(eveny,:,[3 4 1 2]);
%%%