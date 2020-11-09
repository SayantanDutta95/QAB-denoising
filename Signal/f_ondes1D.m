% Sample code of the paper:
% 
% Sayantan Dutta, Adrian Basarab, Bertrand Georgeot, and Denis Kouamé.
% "Quantum mechanics-based signal and image denoising." 
% arXiv preprint arXiv:2004.01078 (2020).
% 
% MATLAB code prepard by Sayantan Dutta
% E-mail: sayantan.dutta@irit.fr and sayantan.dutta110@gmail.com
% 
% This script shows an example of calling our denoising algorithm 
% for signal using Quantum adaptative basis (QAB)

function [psi,E] = f_ondes1D(signal,poids)
% 
%the Hamiltonian associates with the signal and the associated eigenvalues

% Data
N = size(signal,2); 

% creat space to store data
psi = zeros(N,N); %eigenvectors
E = zeros(N,1);   %eigenvalues
        
% Construction of Hamiltonian matrice H
terme_hsm = ones(1,N) * poids;
H = diag(signal,0) + diag(terme_hsm,0)*2 ...
    - diag(terme_hsm(1:N-1),-1) - diag(terme_hsm(1:N-1),1);
H(1,N) = -poids;
H(N,1) = -poids;
    
% Calculation of eigenvalues and eigenvectors
[vectP,valP] = eig(H);
valP = diag(valP);

vp_min = min(valP);
vp_max = max(valP);
    
for g = 1:N
  % Each iteration finds the "following" eigenvector
  %(sorts the vectors in ascending order of the associated eigenvalues)
          
    [valP_assoc,i_psi] = min(valP);
    psi(:,g) = vectP(:,i_psi);
	E(g) = valP(i_psi);
	valP = [valP(1:(i_psi - 1)) ; valP((i_psi + 1):(N - g + 1))];
	vectP = [vectP(:,1:(i_psi - 1)) vectP(:,(i_psi + 1):(N - g + 1))];
        
end
end