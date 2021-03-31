% The annihilating filter method
clear; close;
% B-splines
degree = [0 1 2 3];
N = 2048; % Ngth of kernels of finite support
period = 64; % sampling period
ITER = log2(period); % number of ITERations
shift = 31; % number of shifts
t = 0: 1/period : (N-1)/period; % time of sampling points

[Phi_T] = bspline(period, max(degree)); % bspline reproduces polynomials 
[kernelSet] = kernel(N, period, shift, Phi_T); % obtain kernel 
[dualKernel] = dual_basis(kernelSet(1, :)); % obtain the dual basis kernel
[dualKernelSet] = kernel(N, period, shift, dualKernel); % obtain dual basis kernel set by shifting
[O, R, K] = reproduce(N, period, shift, max(degree), t, dualKernelSet);

% plot the results
figure;
for i = degree
    subplot(2, 2, i+1);
    plot(t, O(i+1, :), 'r','linewidth', 2);
    hold on;
    plot(t, R(i+1, :), 'b', 'linewidth', 2);
    hold on;
    plot(t, K{i+1}, 'y--')
    title(sprintf("Degree: %d", i));
    xlabel('Time');
    ylabel('Amplitude');
end


