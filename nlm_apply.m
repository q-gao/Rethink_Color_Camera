function img = nlm_apply(imin,S,iters,weights)

[h,w,tmp] = size(imin);

Lf = sum(imin,3); imref = imin;

for it = 1:iters
  img = zeros([h w 3]); wts = zeros(h,w);
  idx = 1;
  for i = 0:S
    for j = -S:S
      if i == 0 && j == 0
	continue;
      end;
    
      fprintf('\r [%d] %d / %d       ',it,idx,length(weights));
      wij = double(weights{idx})/255; idx = idx+1;
          
      madd(wij,imref,i,j);
 
      if 1 > 2%
      isp = [1:h-i];
      if j > 0
	jsp = [1:w-j];
      else
	jsp = [(1-j):w];
      end;

      wts(isp+i,jsp+j) = wts(isp+i,jsp+j)+wij;
      wts(isp,jsp) = wts(isp,jsp)+wij;
 
      img(isp,jsp,:) = img(isp,jsp,:) + ...
	  wij(:,:,[1 1 1]).*imref(isp+i,jsp+j,:);
       
      img(isp+i,jsp+j,:) = img(isp+i,jsp+j,:) + ...
	  wij(:,:,[1 1 1]).*imref(isp,jsp,:);
      end;%
      
      
    end;
  end;
  img = img + imref; 
  %wts = wts + 1;  img = img ./ wts(:,:,[1 1 1]);
  Lfz = Lf ./ (sum(img,3)+eps);
  img = img .* Lfz(:,:,[1 1 1]);
  
  imref = img;
end;
fprintf('\n');
