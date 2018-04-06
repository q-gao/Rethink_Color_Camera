% Compute an Edge Map based on luminance image, and normalize per
% 'block' of confusion.
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function emap = EdgeMap(Lf,xidx,yidx,scales)

emap = {0,0,0,0};

Lf = Lf.^0.5;
for i = 1:length(scales)
  
  sgma = scales(i);
  % Use Derivative of Gaussian Filter
  lm = round(4*sgma); x = [-lm:lm]/sgma;
  sx = exp(-x.^2/2);
  fx = sx .* -x/sgma^3;
  
  gx = conv2(conv2(Lf,sx','same'),fx,'same');
  gy = conv2(conv2(Lf,sx,'same'),fx','same');
    
  gxy = (gx+gy).^2/2;
  gyx = (gx-gy).^2/2;
  gy = gy.^2; gx = gx.^2;

  % Thin
  gx = max(0,conv2(gx,[-0.25 1 -0.25],'same'));
  gy = max(0,conv2(gy,[-0.25 1 -0.25]','same'));
  dg = diag([-0.25 1 -0.25]);
  gxy = max(0,conv2(gxy,dg,'same'));
  dg = dg(end:-1:1,:);
  gyx = max(0,conv2(gyx,dg,'same'));
    
  
  % Add each scale in a weighted manner
  emap{1} = emap{1}+gx;
  emap{2} = emap{2}+gy;
  emap{3} = emap{3}+gxy;
  emap{4} = emap{4}+gyx;
end;

for i = 1:4
  emap{i} = emap_norm(emap{i},xidx,yidx);
end;



function emap = emap_norm(em_in,xidx,yidx)

emap = em_in; emean = emap;
for i = 1:max(yidx)
  sub = find(yidx == i);
  emn = mean(emap(sub,:),1);
  emean(sub,:) = repmat(emn,[length(sub) 1]);
end;
for i = 1:max(xidx)
  sub = find(xidx == i);
  emn = mean(emean(:,sub),2);
  emean(:,sub) = repmat(emn,[1 length(sub)]);
end;
emap = emap ./ (eps+emean);
emap = exp(-emap);
