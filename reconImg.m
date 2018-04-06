% reconImg: Reconstruct Image from pattern sampled 
%           by CFA.m
%
% [img,aout,emap] = reconImg(L,f,nzvar)
%
% L: Pattern (from CFA.m)
% f: Color sampling frequency (same as passed to CFA.m)
% nzvar (Optional): Noise Variance
%
% Output:
%
% img:  Reconstructed image
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function img = reconImg(L,f,nzvar)

if nargin == 2
  nzvar = 1e-6;
else
  nzvar = max(nzvar,1e-6);
end;

%%%
fprintf('Getting luminance\n');
Lf = getLum(L,f);

%%%
fprintf('Getting chromaticities\n');
[R,G,B] = getChrom(L,Lf,f,nzvar);

%%%% Fill in missing chromaticities %%%
[xidx,yidx,lxidx,lyidx,rxidx,ryidx] = getIndices(f,size(Lf));

%%%
fprintf('Computing edge map\n');
emap = EdgeMap(Lf,xidx,yidx,[1 2 4]);

%%%
fprintf('Assigning material affinity\n');
alpha = Msegment(emap,xidx,yidx,lxidx,lyidx,rxidx,ryidx);

%%%
fprintf('Combining colors\n');
img = combineColors(alpha,Lf,R,G,B, ...
		    xidx,yidx,lxidx,lyidx,rxidx,ryidx);

img = img * 3;
img = max(0,min(1,img));

img1 = img;
%%%
fprintf('Using non-local matching to improve estimate\n');
img = nlm_chrom(img,f,nzvar);
img = max(0,min(1,img));

