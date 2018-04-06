% medFilt: Median Filter utility function
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>

function [R2,G2,B2] = medFilt(R,G,B)

in = zeros([size(R) 3]); 
in(:,:,1) = R;
in(:,:,2) = G;
in(:,:,3) = B;

% Directional median filters
NFILT = 4; TAP = 5;
odfilt = cell(NFILT,1);
odfilt{1} = ones(1,TAP);
odfilt{2} = ones(TAP,1);
odfilt{3} = diag(ones(TAP,1));
odfilt{4} = odfilt{3}(end:-1:1,:);

outs = zeros([size(in) NFILT]);
TP2 = ceil(TAP/2);

for i = 1:3
  for j = 1:NFILT
    outs(:,:,i,j) = ordfilt2(in(:,:,i),TP2,odfilt{j});
  end;
end;
outs = outs ./ repmat(sum(outs,3),[1 1 3]);

wts = sum(abs(repmat(in,[1 1 1 4])-outs).^2,3);
[tmp,idx] = min(wts,[],4);

out = outs(:,:,:,1);
for j = 2:NFILT
  for i = 1:3
    ci = out(:,:,i);
    co = outs(:,:,i,j);
    ci(idx == j) = co(idx == j);
    out(:,:,i) = ci;
  end;
end;

R2=R; R2(TP2:end-TP2+1,TP2:end-TP2+1) = ...
   out(TP2:end-TP2+1,TP2:end-TP2+1,1);

G2=G; G2(TP2:end-TP2+1,TP2:end-TP2+1) = ...
   out(TP2:end-TP2+1,TP2:end-TP2+1,2);

B2=B; B2(TP2:end-TP2+1,TP2:end-TP2+1) = ...
   out(TP2:end-TP2+1,TP2:end-TP2+1,3);
