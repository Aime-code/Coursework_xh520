% Reconstruction in the presence of noise
clear; close all;
K = 2; % number of diracs
moments = [4, 7, 9]; % number of moments
MAXdegree = moments - 1; % max degree of polynomials
N = 2048; % kernels of finite support
T = 64; % sAling T
maxA = 32; % max Amplitude
shift = 31; % number of shifts
ITER = log2(T); % number of iterations
t = 0: 1/T : (N-1)/T; % time of sAling points
variance = [0.01, 1, 10]; % standard deviation
[signal, location, A] = diracs(N, T, K, maxA); % generate dirac signal

for i_variance = 1:length(variance)
    figure(i_variance);
    noiseSet = sqrt(variance(i_variance)) * randn(1, max(MAXdegree) + 1); % AWGN noise
    fprintf("----Variance = %.2f---- \n", variance(i_variance));
    for m = 1: length(moments) % vary moment
        [phiT, ~, ~] = wavefun(sprintf('dB%d', moments(m)), ITER);
        Kernel = kernel(N, T, shift, phiT);
        [~,~,~,coefs] = reproduce(N, T, shift, MAXdegree(m), t, Kernel);
        samples = signal * Kernel'; % samples
        tau = zeros(1, MAXdegree(m) + 1);
        for i = 0: MAXdegree(m)
            tau(1, i + 1) = dot(coefs(i + 1, :), samples);
        end
        noise = noiseSet(1: MAXdegree(m) + 1);
        Noisy_tau = tau + noise;
        Noisy_TAU = zeros(moments(m) - K, K + 1); % noisy tau matrix
        for j = 1: moments(m) - K
            Noisy_TAU(j, :) = flip(Noisy_tau(j: j + K));
        end

        % Total least squares
        [uTLS, sTLS, vTLS] = svd(Noisy_TAU); % svd decomposition
        hkTLS = vTLS(:, end); 
        location_TLS = sort(real(zero(tf(hkTLS', 1))))'; % determine dirac locations by annihilating filter
        locMatrix = fliplr(vander(location_TLS))'; % Vandermonde system
        TAU3 = tau(1: K)'; 
        A_TLS = (locMatrix \ TAU3)'; % amplitude

        % Cadzow
        [uCadzow, sCadzow, vCadzow] = svd(Noisy_TAU); % svd decomposition
        check = rank(Noisy_TAU) == K; % tau matrix is full rank in the noisy case
        toeplitz = zeros(size(Noisy_TAU));
    
        % Iterate tau matrix 
        while ~check
            sCadzow(size(sCadzow, 2), size(sCadzow, 2)) = 0;
            % reconstruct tau matrix
            Noisy_TAU = uCadzow * sCadzow * vCadzow';
            % it is not Toeplitz now; average diagonals to make it Toeplitz
            for i = 1: size(Noisy_TAU, 1)
                for j = 1: size(Noisy_TAU, 2)
                toeplitz(i, j) = mean(diag(Noisy_TAU, j - i));
                end
            end
            Noisy_TAU = toeplitz;
            [uCadzow, sCadzow, vCadzow] = svd(Noisy_TAU);
            check = rank(sCadzow) == K;
        end

        hkCadzow = vCadzow(:, end);
        location_Cadzow = sort(real(zero(tf(hkCadzow', 1))))';
        locMatrix = fliplr(vander(location_Cadzow))'; % Vandermonde system
        TAU3 = tau(1: K)'; 
        A_Cadzow = (locMatrix \ TAU3)'; % amplitude

        % Plot the results
        sgtitle(sprintf('variance = %.2f', variance(i_variance)));
        subplot(3,1,m);
        stem(location, A, 'k');
        hold on;
        stem(location_TLS, A_TLS, 'r', 'linewidth', 1.5);
        hold on;
        stem(location_Cadzow, A_Cadzow, 'b');
        hold on;
        xlabel('Time');
        ylabel('Amplitude');
        legend('Original', 'TLS', 'Cadzow');
        title(sprintf('Moments = %d', moments(m)));
        fprintf("Moment = %d \n", moments(m));
        
        fprintf("Real location: %.5f and %.5f \n", location);
        fprintf("Real Amplitude: %.5f and %.5f \n", A);
        fprintf("Estimated TLS location: %.5f and %.5f \n", location_TLS);
        fprintf("Estimated TLS Amplitude: %.5f and %.5f \n", A_TLS);
        fprintf("Estimated Cadzow location: %.5f and %.5f \n", location_Cadzow);
        fprintf("Estimated Cadzow Amplitude: %.5f and %.5f \n\n", A_Cadzow);
    end
end
   