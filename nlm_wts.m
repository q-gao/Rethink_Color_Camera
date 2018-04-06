function weights = nlm_wts(imin,S,K,sgm)

% S = Search window size = 2S+1
% K = Match window size = 2K+1

weights = cell((S+1)*(2*S+1)-1,1);

[h,w,tmp] = size(imin);
sgm = 2*sgm.^2;

kflt = fspecial('gaussian',[(2*K+1) 1],K/2);
kflt = kflt / sqrt(sum(kflt(:).^2));

Lf = sum(imin,3);

idx = 1;
for i = 0:S
  for j = -S:S
    if i == 0 && j == 0
      continue;
    end;
    
    fprintf('\r %d / %d       ',idx,length(weights));
    
    isp = [1:h-i];
    if j > 0
      jsp = [1:w-j];
    else
      jsp = [(1-j):w];
    end;
    
    wij = (Lf(isp,jsp)-Lf(isp+i,jsp+j)).^2;
    wij = conv2(conv2(wij,kflt','same'),kflt,'same');
    wij = exp(-wij./sgm);
    
    weights{idx} = uint8(wij*255);
    idx = idx+1;
  end;
end;
fprintf('\n');