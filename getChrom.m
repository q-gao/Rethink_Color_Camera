% getChrom: Get chromaticity utility function
%
% Assumes same "material" in each color sample block. Also returns an
% initial "smoothed" estimate of luminance, that can be used to for
% boundary blocks ignored by getLum.
%
% -- Ayan Chakrabarti <ayanc@eecs.harvard.edu>
function [R,G,B,Lf] = getChrom(L,Lf,f,nzvar)

X11 = zeros(size(L(1:f:end,1:f:end)));
X12 = X11; X22 = X11;
Y1 = X11; Y2 = X22;

sz = size(L);
x1 = [1:f:sz(2)]; x2 = [2:f:sz(2)];
y1 = [1:f:sz(1)]; y2 = [2:f:sz(1)];
lx1 = [1:length(x1)]; lx2 = [1:length(x2)];
ly1 = [1:length(y1)]; ly2 = [1:length(y2)];

% Regularizer
X11(:) = nzvar; X22 = X11;
Y1(:) = nzvar/3; Y2 = Y1;

% R 
X11(ly1,lx2) = X11(ly1,lx2) + Lf(y1,x2).^2;
Y1(ly1,lx2) = Y1(ly1,lx2) + Lf(y1,x2).*L(y1,x2);

% G1
X22(ly1,lx1) = X22(ly1,lx1) + Lf(y1,x1).^2;
Y2(ly1,lx1) = Y2(ly1,lx1) + Lf(y1,x1).*L(y1,x1);

% G2
X22(ly2,lx2) = X22(ly2,lx2) + Lf(y2,x2).^2;
Y2(ly2,lx2) = Y2(ly2,lx2) + Lf(y2,x2).*L(y2,x2);

% B
l = Lf(y2,x1).^2; r = Lf(y2,x1).*max(0,Lf(y2,x1) - L(y2,x1));
X11(ly2,lx1) = X11(ly2,lx1) + l;
X22(ly2,lx1) = X22(ly2,lx1) + l;
X12(ly2,lx1) = X12(ly2,lx1) + l;

Y1(ly2,lx1) = Y1(ly2,lx1) + r;
Y2(ly2,lx1) = Y2(ly2,lx1) + r;

det = 1./(X11.*X22 - X12.^2);
R = det.*(X22.*Y1 - X12.*Y2);
G = det.*(X11.*Y2 - X12.*Y1);
B = max(0,1-R-G);

% Normalize just in case
R = min(1,R); G = min(1,G); B = min(1,B);
sm = R+G+B;
R = R ./ sm;G = G ./ sm;B = B ./ sm;
R(isnan(R)) = 0; G(isnan(G)) = 0; B(isnan(B)) = 0;

for it=1:50
  [R,G,B] = medFilt(R,G,B);
end;
