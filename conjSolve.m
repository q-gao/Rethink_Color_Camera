function a = conjSolve2(a0,emap,msk,numiter,tol);

a = a0;
for i = 1:size(a0,3)
  a(:,:,i) = csolve(a0(:,:,i),emap,msk,numiter,tol);
end;

function a = csolve(a0,emap,msk,numiter,tol)

wx = min(emap{1}(:,1:end-1),emap{1}(:,2:end));
wy = min(emap{2}(1:end-1,:),emap{2}(2:end,:));
wxy = min(emap{3}(1:end-1,1:end-1),emap{3}(2:end,2:end));
wyx = min(emap{4}(1:end-1,2:end),emap{4}(2:end,1:end-1));

r = Apx(a0,wx,wy,wxy,wyx,msk);
r = -r;
p = r;
a = a0;

adiff = zeros(size(a));
for iters=1:numiter
  fprintf('\r[%02d] Cost = %e         ',iters,sum(r(:).^2));
  Ap = Apx(p,wx,wy,wxy,wyx,msk);
  alpha = sum(r(:).^2) / sum(p(:).*Ap(:));
  a = a + alpha * p; rnew = r-alpha*Ap;
  if mean(r(:).^2) < tol
    break;
  end;
  beta = sum(rnew(:).^2)/sum(r(:).^2);
  r = rnew; p = r+beta*p;
end;
fprintf('\n');
a = max(0,min(1,a));

function Ap = Apx(p,wx,wy,wxy,wyx,msk)
gf = [-1 1]; gd1 = diag(gf); gd2 = gd1(end:-1:1,:);

Ap=conv2(conv2(p,gf,'valid').*wx,-gf,'full');
Ap=Ap+conv2(conv2(p,gf','valid').*wy,-gf','full');
Ap=Ap+conv2(conv2(p,gd1','valid').*wxy,-gd1','full');
Ap=Ap+conv2(conv2(p,gd2','valid').*wyx,-gd2','full');
Ap = Ap.*msk;