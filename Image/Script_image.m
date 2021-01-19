% Sample code of the paper:
% 
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé.
% "Quantum mechanics-based signal and image representation: application to denoising." 
% arXiv preprint arXiv:2004.01078 (2020).
% 
% MATLAB code prepard by Sayantan Dutta
% E-mail: sayantan.dutta@irit.fr and sayantan.dutta110@gmail.com
% 
% This script shows an example of calling our denoising algorithm 
% for image using Quantum adaptative basis (QAB)

close all
clc; clear;

%load data
load('sample_image.mat')

%size of signal
[M,N] = size(I);
NN = N^2;
%generate the noise
SNR = 15;
pI = sum(I(:) .^2) / length(I(:));
B = randn(N,N) .* sqrt(abs(I));  % Poisson noise
pB_tmp = sum(sum(B .^2)) / NN;
B = B / sqrt(pB_tmp) * sqrt(pI * 10^(- SNR / 10));
J = B + I;  %noisy image

% Data
Ms = 20; % number of iteration in reconstruction
pds = 3; % value of palnck's constant
sg = 7.5; % Gaussian Variance (smoothing)

%image denoising using QAB
[I_result] = image_denoising_QAB(I,J,Ms,pds,sg);

% Here one may think why we need I as input, since I is the clean image.
% Actually in the program, the clean image data I is used just to find out the 
% best thresholding hyperparameter. In the code, one can see that the clean 
% image I is only used for the computation of signal-to-noise-ratio (SNR) 
% and based on this one can tune the best thresholding hyperparameter.
% So, one does not need any knowledge about the clean image I for computing
% the quantum adaptive basis (QAB) and can tune the thresholding hyperparameter
% manually up to their best possible values and modify the code accordingly.
% Therefore, one should not get confused after seeing I as an input, because I
% does not have any impact on the denoising process.


% Calculation of SNR, PSNR and SSIM
pnB = sum((I_result(:) - I(:)) .^2) / length(I(:));
SNR_end = 10 * log10(pI / pnB);
PSNR_end = calc_PSNR(I,I_result);
SSIM_end = ssim_index(I_result, I, [0.01 0.03], fspecial('gaussian', 3, 1.5), max(I(:)))
% Print data
fprintf('\n\n OUTPUT:\n SNR = %.2f, PSNR = %.2f, and SSIM = %.2f \n\n',SNR_end,PSNR_end,SSIM_end);

%plot figures
font_size = 12;
figure;
subplot(1,3,1);
imagesc(I); colormap gray
xlabel('Horizontal distance [pixels]','Fontsize',font_size)
ylabel('Horizontal distance [pixels]','Fontsize',font_size)
set(gca,'FontSize',font_size);
title(sprintf('Clean image'))

subplot(1,3,2);
imagesc(J); colormap gray
xlabel('Horizontal distance [pixels]','Fontsize',font_size)
ylabel('Horizontal distance [pixels]','Fontsize',font_size)
set(gca,'FontSize',font_size);
title(sprintf('Noisy image (SNR = %.2f dB)',SNR))

subplot(1,3,3);
imagesc(I_result); colormap gray
xlabel('Horizontal distance [pixels]','Fontsize',font_size)
ylabel('Horizontal distance [pixels]','Fontsize',font_size)
set(gca,'FontSize',font_size);
title(sprintf('Denoised image(PSNR = %.2f dB)',PSNR_end))
