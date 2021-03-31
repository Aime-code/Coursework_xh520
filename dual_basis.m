function [dualKernel] = dual_basis(kernel)

% obtain the dual basis kernel of a biorthogonal wavelet
% Input:
% kernel: basis kernel of a biorthogonal wavelet
% Output:
% dual: the dual basis kernel of a biorthogonal wavelet

gram = zeros(length(kernel)); % Gram matrix
for i = 1: length(kernel)
    for j = 1: length(kernel)
       gram(i, j) =  dot(kernel(i), kernel(j));
    end
end
gram = gram / norm(gram); 
dualKernel = kernel * gram; % dual basis
end

