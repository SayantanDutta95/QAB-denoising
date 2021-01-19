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
% for signal using Quantum adaptative basis (QAB)

close all
clear;clc;

%load data
load('sample_signal.mat');

%size of signal
N = size(S,2); 
%generate the noise
SNR = 15;
pS = sum(S .^2) / N;
B = randn(1,N) .* sqrt(abs(S)); % Poisson noise
pB_tmp = sum(sum(B .^2)) / N;
B = B / sqrt(pB_tmp) * sqrt(pS * 10^(- SNR / 10));
SB = B + S; %noisy signal

% data
Ms = 200; % number of iteration in reconstruction
pds = 0.4; % value of palnck's constant
sg = 15; % Gaussian Variance (smoothing)

%signal denoising using QAB
[S_result] = signal_denoising_QAB(S,SB,Ms,pds,sg);

% Here one may think why we need S as input, since S is the clean signal.
% Actually in the program, the clean signal data S is used just to find out the 
% best thresholding hyperparameter. In the code, one can see that the clean 
% signal S is only used for the computation of signal-to-noise-ratio (SNR) 
% and based on this one can tune the best thresholding hyperparameter.
% So, one does not need any knowledge about the clean signal S for computing
% the quantum adaptive basis (QAB) and can tune the thresholding hyperparameter
% manually up to their best possible values and modify the code accordingly.
% Therefore, one should not get confused after seeing S as an input, because S
% does not have any impact on the denoising process.


%Calculation of SNR and PSNR
pnB = sum((S_result(:) - S(:)) .^2) / length(S(:));
SNR_end = 10 * log10(pS / pnB);
PSNR_end = calc_PSNR(S,S_result);
% Print data
fprintf('\n\n OUTPUT:\n SNR = %.2f and PSNR = %.2f \n\n',SNR_end,PSNR_end);

%plot figure
figure;
subplot(1,3,1);plot(S);axis([0 512 -10 25]);
title(sprintf('Clean signal'))
subplot(1,3,2);plot(SB);axis([0 512 -10 25]);
title(sprintf('Noisy signal (SNR = %.2f dB)',SNR))
subplot(1,3,3);plot(S_result);axis([0 512 -10 25]);
title(sprintf('Denoised signal (PSNR = %.2f dB)',PSNR_end));
    
