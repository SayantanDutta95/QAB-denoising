function H = heavi(x,eps)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
N = size(x,2);
fonct = zeros(1,N);

gauche = ones(1,N);
milieu = 0.5 * (1 - x / eps - 1 / pi * sin(pi * x / eps));

H = fonct + gauche .* (x <= -eps) + milieu .* (x > -eps) .* (x < eps);

end

