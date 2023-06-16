%% load image
load kidney.mat

%% take user input and calculate geodesic distance
[gd,mask] = geodist(im);

%% run first stage optimisation
lambda = 5;
mu = 1;
eta = 20; %tweak this parameter

u = SelectiveMS_PrimalDual(im,lambda,mu,eta,gd);

figure;
subplot(1,2,1); imagesc(im); axis image; colormap gray; title("Input image");
subplot(1,2,2); imagesc(u); axis image; colormap gray; title("Model one output");
%% segment via thresholding
th=0.35; %tweak this parameter
seg = zeros(size(u));  seg(u>th)=1;

figure; imagesc(seg); colormap gray; axis image; title("Binary output");
figure; imagesc(im); colormap gray; axis image; hold on; 
contour(seg,[0.5,0.5],'r','LineWidth',2); title("Contour plot");

