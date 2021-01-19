function PSNR = calc_PSNR(img1,img2)
%% US_ADM_calc_PSNR Computes the Peak Signal-to-Noise Ratio.
%  PSNR = calc_PSNR(img1,img2)
%
%   img1: reference image
%   img2: output image
%
%  ----------------------
% | Code:  Renaud Morin  |
% | Email: morin@irit.fr |
% | Date:  July 2011     |
%  ----------------------


%% Initialization
if nargin~=2
    error('MATLAB:paramAmbiguous','There must be exactly 2 arguments.')
end
if size(img1)~=size(img2)
    error('MATLAB:paramAmbiguous','Inputs must be of same size')
end


%% Process 
d= max([img1(:);img2(:)]);
mse=sum( (img1(:)-img2(:)).^2 )/numel(img1);

PSNR=10*log10(d^2/mse);
