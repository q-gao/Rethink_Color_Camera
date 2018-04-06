% getLum: Get Luminance utility function
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function Lf = getLum(L,f,nzvar)

% Mask
lmask = ones(size(L));
lmask(1:f:end,2:f:end) = 0; lmask(2:f:end,1:f:end) = 0;
lmask(1:f:end,1:f:end) = 0; lmask(2:f:end,2:f:end) = 0;

betas = 2.^(-0.25*[0:49]);

Lf = L;
for i = 1:length(betas) 
  fprintf('\r   Shrinkage step %d of %d       ',i,length(betas));
  Lf = shrinkL(Lf,betas(i));
  Lf = (1-lmask).*Lf + lmask.*L;
end;
fprintf('\n');
  
Lf = max(0,Lf);

% Do a wavelet decomposition and shrink
function Lf = shrinkL(L,b)

%h = [1 1]/sqrt(2); g = [1 -1]/sqrt(2);
[h,g] = wfilters('db2');

% Decompose
coeffs = cell(3,1); c = L;
gx = conv2(c,g(end:-1:1),'full');
hx = conv2(c,h(end:-1:1),'full');
  
% gxhy
coeffs{1} = conv2(gx,h(end:-1:1)','full');
% hxgy
coeffs{2} = conv2(hx,g(end:-1:1)','full');
% gxgy
coeffs{3} = conv2(gx,g(end:-1:1)','full');
% Scaling
Lf = conv2(hx,h(end:-1:1)','full');

% Shrink
for j = 1:3
  tv2 = sqrt(conv2(coeffs{j}.^2,ones(3,3)/9,'same'));
  scale = max(0,tv2-b) ./ (tv2 + eps);
  coeffs{j} = coeffs{j} .* scale;
end;


% Reconstruct
Lf = spConv(Lf,h'*h,6);
Lf = Lf + spConv(coeffs{1},h'*g,6);
Lf = Lf + spConv(coeffs{2},g'*h,6);
Lf = Lf + spConv(coeffs{3},g'*g,6);

Lf = Lf /4;
