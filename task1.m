% Strang-Fix conditions
clear; close;
degree = [0 1 2 3];
N = 2048; % the N of kernels of finite support
period = 64; % the sampling period
shift = 31; % the number of shifts
ITER = log2(period); % the number of iterations
t = 0: 1/period : (N-1)/period; % time of sampling points

[phi_T, ~, ~] = wavefun('dB4', ITER);
[shift_phi] = kernel(N, period, shift, phi_T); % obtain kernel by shifting scaling function
[O, R, K] = reproduce(N, period, shift, max(degree), t, shift_phi); % determine polynomials and coefficients of corresponding kernels

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