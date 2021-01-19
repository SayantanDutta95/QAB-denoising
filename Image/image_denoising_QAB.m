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

function [I_result] = image_denoising_QAB(I,J,Ms,pds,sg)

% Here one may think why we need I as input, since I is the clean image.
% Actually in the program, the clean image data I is used just to find out the 
% best thresholding hyperparameter. In the code, one can see that the clean 
% image I is only used just for the computation of signal-to-noise-ratio (SNR) 
% and based on this one can tune the best thresholding hyperparameter.
% So, one does not need any knowledge about the clean image I for computing
% the quantum adaptive basis (QAB) and can tune the thresholding hyperparameter
% manually up to their best possible values and modify the code accordingly.
% Therefore, one should not get confused after seeing I as an input, because I
% does not have any impact on the denoising process.


%size of image
[M,M1] = size(I); % assume M = M2
MM = M^2;
N = 64;
NN = N^2;
P = 2*M/N;
% store initial data
I_old = I;
J_old = J;

saut = 12; %Moving the threshold (must be at least 1, this is the slowest case)

%  Threshold parameter
Vs = linspace(7,11,Ms); % for reconstruction
Vs = 2 .^Vs;      

disp('Start the search of wave function')
J_new = zeros(M,M);  %creat space for the big image
o = 0;
cmpt = zeros(M,M);   %count the overlapping

for i = 0:(P-2)
    for j = 0:(P-2)
        
   % These variables store the values associated with the maximums.
    V_ms = 0.01;
    RSB_ms = 0; % maximum SNR for adaptive transformation 
    seuil_ms = 1;
        
  % Devide the image in to small parts
        J_part = J((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N));
        I_part = I((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N));

  % Normalize small part
        J_part_max = max(max(J_part));
        I_part = I_part/J_part_max;
        J_part = J_part/J_part_max;
   
        pI_part = sum(sum(I_part .^2)) / NN;

%  The noisy image is blurred
%  Creation of a Gaussian
[x,y] = meshgrid((-N/2):(N/2 - 1));
z = 1 / (sqrt(2 * pi * sg))^2 * exp(-(x .^2 + y .^2) / (2 * sg));

% Convolution product
gaussF = fft2(ifftshift(z));
JF = fft2(J_part);
n_J = real(ifft2(gaussF .* JF));

% new images
Jns = zeros(N,N);
%Calculation of the eigenvectors then projection in the base
J_col = reshape(J_part,NN,1);
[psi,psi_col,E] = f_ondes2D(n_J,pds);
 alp = linsolve(psi_col,J_col);

% Reconstruction des images
for k = 1:Ms %Loop on the slope parameter
    
    v = Vs(k);
        
    max_atteint = 0;
    RSB_s = 0;
    l = 1;
        
    while (l <= NN) * (max_atteint < 20) %Loop on the threshold position, stops when one reaches a local maximum (global?) For the threshold parameter
      
        % threshold
        x = (1:(NN)) - l + 2;
        taux = heavi(x,v);
        
        %Re-construction
        n_I = zeros(N,N);
        for t = 1:(NN)
            n_I = n_I + psi(:,:,t) * taux(t) * alp(t);
        end
    
        %  Calculation of SNR
        n_B = n_I - I_part;
        pnB = sum(sum(n_B .^2)) / NN;
        RSB_n = 10 * log10(pI_part / pnB);
        
        % Test if the SNR increases
        max_atteint = (max_atteint + 1) * (RSB_s > RSB_n);
        RSB_s = RSB_s * (RSB_s > RSB_n) + RSB_n * (RSB_s <= RSB_n);
        
        V_ms = V_ms * (RSB_ms >= RSB_n) + v * (RSB_ms < RSB_n);
        seuil_ms = seuil_ms * (RSB_ms >= RSB_n) + l * (RSB_ms < RSB_n);
        Jns = Jns * (RSB_ms >= RSB_n) + n_I * (RSB_ms < RSB_n);
        RSB_ms = RSB_ms * (RSB_ms >= RSB_n) + RSB_n * (RSB_ms < RSB_n);
                
        l = l + saut;
    end
end

%store the small part into the big image
        J_new((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N)) = Jns * J_part_max + J_new((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N));
        cmpt((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N)) = cmpt((i*N/2+1):(i*N/2+N),(j*N/2+1):(j*N/2+N)) + 1; 
        o = o+1
        
    end
end

J_new = J_new./cmpt; %remove the overlapping

I_result = J_new;

end