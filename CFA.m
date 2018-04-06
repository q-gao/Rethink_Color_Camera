% CFA: Make a sampled single channel image, based on
%      sparse color samples. Note that the colors are
%      sampled as per the regular Bayer [GR;BG] pattern.
%
% [L,Limg] = CFA(img,f)
% 
% img: RGB Image
% f:   Color sampling interval (must be greater than 2)
%
% Output:
% L:   Sampled pattern (scaled down by 3 to maintain range)
% Limg: 3-color Visualization of L
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function [L,Limg] = CFA(img,f)

L = sum(img,3);
L(1:f:end,1:f:end) = img(1:f:end,1:f:end,2);
L(2:f:end,2:f:end) = img(2:f:end,2:f:end,2);

L(1:f:end,2:f:end) = img(1:f:end,2:f:end,1);
L(2:f:end,1:f:end) = img(2:f:end,1:f:end,3);

L = L / 3;

if nargout == 2
  Limg = visCFA(L,f);
end;

function img = visCFA(L,f)

img = repmat(L,[1 1 3])*3;

img(1:f:end,1:f:end,[1 3]) = 0;
img(2:f:end,2:f:end,[1 3]) = 0;
img(1:f:end,2:f:end,[2 3]) = 0;
img(2:f:end,1:f:end,[1 2]) = 0;
