% Form full image based on colors at corners and material
% affinities.
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function img = combineColors(alpha,Lf,R,G,B,xidx,yidx,lxidx,lyidx,rxidx,ryidx)

% Get luminance at four corners
Lij = zeros([size(Lf) 4]);
Lij = Lf(lyidx,lxidx); Lij(:,:,2) = Lf(lyidx,rxidx);
Lij(:,:,3) = Lf(ryidx,lxidx); Lij(:,:,4) = Lf(ryidx,rxidx);

% Go from alpha -> beta by matching luminance
beta = alpha.^2;
nrm = sum(beta.*Lij,3) ./ (Lf+eps);
beta = beta ./ repmat(nrm+eps,[1 1 4]);

% Chromaticities at four corners
R = [R R(:,end)]; R = [R; R(end,:)];
G = [G G(:,end)]; G = [G; G(end,:)];
B = [B B(:,end)]; B = [B; B(end,:)];

% Do this so that beta*R is the contribution of Red intensity and
% not Red chromaticity.
beta = beta .* Lij;

% Combine
img = zeros([size(Lf) 3]);
k=1;
for dy = 0:1
  for dx = 0:1
    img(:,:,1)=img(:,:,1) + beta(:,:,k).*R(yidx+dy,xidx+dx);
    img(:,:,2)=img(:,:,2) + beta(:,:,k).*G(yidx+dy,xidx+dx);
    img(:,:,3)=img(:,:,3) + beta(:,:,k).*B(yidx+dy,xidx+dx);
    
    k = k+1;
  end;
end;
