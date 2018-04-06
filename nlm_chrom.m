function imr = nlm_chrom(img,f,nzvar)

if nargin == 2
  nzvar = 1e-6;
end;
nzvar = max(1e-6,nzvar);

fprintf('Computing weights\n');
weights = nlm_wts(img,ceil(f/2),ceil(f/4),40*sqrt(nzvar));
fprintf('Applying non-local kernel\n');
imr = nlm_apply(img,ceil(f/2),5,weights);

imr = max(0,min(1,imr));