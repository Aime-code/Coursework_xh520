function [O, R, K, coefficients] = reproduce(N, period, shift, maxdegree, t, shift_phi)
% Function: 
% create polynomials of order up to a maximum number
% calculate the coefficients that the shifted kernels need to reconstruct the corresponding polynomials
% Input:
% N: signal length
% period: sampling period
% shift: number of available shifts in the signal range
% maxdegree: maximum degree of the polynomialsnomials
% t: sampling points
% shift_phi: shifted kernel sets for sampling
% Output:
% O: original polynomials 
% R: reproduced polynomials 
% K: shifted kernels

O = zeros(maxdegree+1, N);
K = cell(maxdegree+1);
coefficients = zeros(maxdegree+1, shift + 1);
for i = 0: maxdegree
    O(i+1, :) = t.^i;     
    coefficients(i+1, :) = dot(repmat(O(i+1, :), shift+1, 1), shift_phi, 2).' ./ period; % coefficients of corresponding kernels
    K{i+1} = shift_phi' .* coefficients(i+1, :);
end
R = coefficients * shift_phi;
end

