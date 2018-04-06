nz_sigma = 2e-4; % Simulated noise std.
f =  6;          % Color sampling rate


% Read ground truth image
img = im2double(imread('test_image.png'));

% Simulate Sampling
L = CFA(img,f);

% Simulate adding noise
L = max(0,min(1,L + randn(size(L)) * nz_sigma));

% Call re-construction algorithm
img_recon = reconImg(L,f,nz_sigma^2);

% Show
imshow(img_recon);

mse = sqrt(mean((img(:)-img_recon(:)).^2));
psnr = 20*log10(1/mse);
fprintf('PSNR = %.2f dB\n',psnr);