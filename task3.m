%  Sampling Diracs
clear; close;
% Annihilating filter
load('tau.mat');
MAXdegree = 3; % max degree of polynomials
K = 2; % number of pulses
TAU1 = zeros(MAXdegree-K + 1, K); % left tau matrix
TAU2 = zeros(MAXdegree-K + 1, 1); % right tau matrix

% Yule-Walker 
for i = 1: MAXdegree-K+1
   TAU1(i, :) = flip(tau(i: i+K-1));
   TAU2(i) = -tau(i + K);
end

% solve equations for filter coefficients
hk = [1; TAU1\TAU2]; % h is called annihilating filter of τ
equation = @(location) ([location(1)*location(2)-hk(3); -(location(1)+location(2))-hk(2)]);
tk = sort(fsolve(equation, [1 1])); % location
tks = fliplr(vander(tk))'; % Vandermonde system
TAU3 = tau(1: K)'; % tau matrix in the Vandermonde system
ak = (tks \ TAU3)'; % amplitude

% Plot Results
stem(tk, ak);
xlabel('Time');
ylabel('Amplitude');
title('Reconstruction of Dirac Signal');
fprintf("Estimated location: %.4f and %.4f \n", tk);
fprintf("Estimated amplitude: %.4f and %.4f \n", ak);
