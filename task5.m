%  Sampling Diracs
clear; close;
load ('samples.mat');
K = 2; % number of diracs
MAXdegree = ceil(2 * K - 1); % max degree of polynomials
N = 2048; % kernels of finite support
T = 64; % sampling T
maxA = 32; % max amplitue
shift = 31; % number of shifts
ITER = log2(T); % number of Iterations
t = 0: 1/T : (N-1)/T; % time of sampling points

[phiT, ~, ~] = wavefun('dB4', ITER);
% Recover signal from samples
tau = zeros(1, MAXdegree + 1);
TAU1 = zeros(MAXdegree-K+1, K); % left tau Matrix
TAU2 = zeros(MAXdegree-K+1, 1); % right tau Matrix
Kernel = kernel(N, T, shift, phiT);
[~,~,~,coefs] = reproduce(N, T, shift, MAXdegree, t, Kernel);

for i = 0: MAXdegree
    tau(1, i + 1) = dot(coefs(i + 1, :), y_sampled);
end
% Yule-Walker 
for j = 1: MAXdegree-K+1
   TAU1(j, :) = flip(tau(j: j + K - 1));
   TAU2(j) = -tau(j + K);
end
hk = [1; TAU1 \ TAU2];
tk = sort(zero(tf(hk',1)))'; % location
tks = fliplr(vander(tk))'; % Vandermonde system
TAU3 = tau(1: K)'; % tau matrix in the Vandermonde system
ak = (tks \ TAU3)'; % Amplitude

% Plot the results
stem(tk, ak, 'o');
xlabel('Time');
ylabel('Amplitude');
title('Reconstruction of Dirac Signal');
fprintf("Estimated location: %.5f and %.5f \n", tk);
fprintf("Estimated amplitude: %.5f and %.5f \n", ak);
