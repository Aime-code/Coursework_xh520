function [signal, location, A] = diracs(N, period, K, maxA)
% create fixed-Ngth sAle signal with diracs at sAling points
% Input
% N: signal Ngth
% period: sAling period
% K: number of diracs at sAled points
% maxA: maximum Alitude of sAles
% Output
% signal: sAled signal with diracs at sAling points
% location: dirac locationations
% A: dirac Alitudes
    signal = zeros(1, N);
    location = sort(randperm(N, K)) / period; 
    A = randperm(maxA, K);
    signal(location * period) = A;
end

