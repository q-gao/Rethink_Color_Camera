function img = nlm_dnz(imin,nzvar,S,K)

if nargin < 3
  S = 8; K = 2;
end;

if nargin == 1
  nzvar = 1e-6;
else
  nzvar = max(nzvar,1e-6);
end;
sgm = 16*nzvar;


[h,w,tmp] = size(imin);

%%% Color Space transform
C = [1;1;1];
C = [C null(C')/sqrt(3)];
imref = reshape(reshape(imin,[h*w 3])*C,[h w 3]);
%%%

sgm = 2*sgm;

kflt = fspecial('gaussian',[(2*K+1) 1],K/2);
kflt = kflt / sqrt(sum(kflt(:).^2));

img = zeros([h w 3]);
wts = zeros([h w]);

idx = 1;
for i = 0:S
  for j = -S:S
    if i == 0 && j == 0
      continue;
    end;
    
    fprintf('\r %d / %d      ',idx,((S+1)*(2*S+1)-1));
    idx = idx + 1;
    
    isp = [1:h-i];
    if j > 0
      jsp = [1:w-j];
    else
      jsp = [(1-j):w];
    end;
    
    im1 = imref(isp,jsp,:);
    im2 = imref(isp+i,jsp+j,:);
    
    wij = sum((im1-im2).^2,3);
    wij = conv2(conv2(wij,kflt','same'),kflt,'same');
    wij = exp(-wij./sgm);

    madd(wij,imref,i,j);

  end;
end;
img = img + imref; wts = wts + 1;
img = img ./ repmat(wts,[1 1 3]);

C = inv(C);
img = reshape(reshape(img,[h*w 3])*C,[h w 3]);
fprintf('\n');