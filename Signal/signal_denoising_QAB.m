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

function [S_result] = signal_denoising_QAB(S,SB,Ms,pds,sg)

% Here one may think why we need S as input, since S is the clean signal.
% Actually in the program, the clean signal data S is used just to find out the 
% best thresholding hyperparameter. In the code, one can see that the clean 
% signal S is only used just for the computation of signal-to-noise-ratio (SNR) 
% and based on this one can tune the best thresholding hyperparameter.
% So, one does not need any knowledge about the clean signal S for computing
% the quantum adaptive basis (QAB) and can tune the thresholding hyperparameter
% manually up to their best possible values and modify the code accordingly.
% Therefore, one should not get confused after seeing S as an input, because S
% does not have any impact on the denoising process.


%normalized the signal
S_i = S; %% store initial data
SB_i = SB;
SB_max = max(abs(SB));
S = S/SB_max;
SB = SB/SB_max;
%size of signal
N = size(S,2); 
saut = 2; % Threshold displacement (must be at least 1, this is the slowest case)
pS = sum(S .^2) / N;

% thresholding parameter we want to make minimizers
Vs = linspace(1,10,Ms); % for reconstruction
Vs = 2 .^Vs;

% These variables store the values associated with the maximums.
V_ms = 0.01;
seuil_m = 1;
RSB_ms = 0; % maximum SNR for adaptive transformation 
S_s_m = zeros(1,N);
pds_m = 10;

disp('simulation using adaptive transformation')

% Gaussian smoothing
x = (-N/2):(N/2 - 1);
x = [x((N/2 + 1):N) x(1:(N/2))];
y = 1 / sqrt(2 * pi * sg) * exp(-x .^2 / (2 * sg ));
SB_f = fft(SB) ;
y_f = fft(y);
SB_flou = real(ifft(SB_f .* y_f));

% Reconstruction of signal
    
    % Calculation of projection using eigen functions
    [psi,E] = f_ondes1D(SB_flou,pds);
    SB_c = transpose(SB);
    a = linsolve(psi,SB_c);
    
    compteur = 1;
    for v = Vs
    	disp(sprintf('Threshold slope parameter: %.0f/%.0f - palnck = %.1f',compteur,Ms,pds));
        for k = 1:2:N
        
        	% Threshold
            x = (1:N) - k + 2;
            taux = heavi(x,v);
        
        	%Re-construction
            n_SB = zeros(1,N);
            for j = 1:N
                n_SB = n_SB + taux(j) * a(j) * transpose(psi(:,j));
            end
        
        	%  Calculation of SNR
        	n_B = n_SB - S;
        	pnB = sum(n_B .^2) / N;
        	RSB_n = 10 * log10(pS / pnB);
    
         	%Test if SNR increases
         	V_ms = V_ms * (RSB_ms >= RSB_n) + v * (RSB_ms < RSB_n);
         	seuil_m = seuil_m * (RSB_ms >= RSB_n) + k * (RSB_ms < RSB_n);
        	S_s_m = S_s_m * (RSB_ms >= RSB_n) + n_SB * (RSB_ms < RSB_n);
            pds_m = pds_m * (RSB_ms >= RSB_n) + pds * (RSB_ms < RSB_n);
            RSB_ms = RSB_ms * (RSB_ms >= RSB_n) + RSB_n * (RSB_ms < RSB_n);
        
        end
        compteur = compteur + 1;
    end
    
S = S_i;                        % original signal
SB = SB_i;                      % noisy signal
S_result = S_s_m.*SB_max;       % denoised signal
end
