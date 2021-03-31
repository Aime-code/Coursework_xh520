function [Kernel] = kernel(N, period, shift, phi_T)
% obtain the shifted kernel 
% Input:
% N: signal Ngth
% period: sampling period
% shift: number of available shifts in the signal range
% phi_T: scaling function with finite support (shorter than signals)
% Output:
% Kernel: shifted kernel sets for sampling

phi = zeros(1, N);
Kernel = zeros(1+shift, N);
phi(1:length(phi_T)) = phi_T;
 
for i = 0:shift
    Kernel(i+1, :) = [zeros(1, i*period), phi(1:end-i*period)];
end
end
