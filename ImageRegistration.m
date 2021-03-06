function [Tx_RGB, Ty_RGB]= ImageRegistration
% *************************************************************************
% Wavelets and Applications Course - Dr. P.L. Dragotti
% MATLAB mini-project 'Sampling Signals with Finite Rate of Innovation'
% Exercice 6
% *************************************************************************
%
% FOR STUDENTS
%
% This function registers the set of 40 low-resolution images
% 'LR_Tiger_xx.tif' and returns the shifts for each image and each layer
% Red, Green and Blue. The shifts are calculated relatively to the first
% image 'LR_Tiger_01.tif'. Each low-resolution image is 64 x64 pixels.
%
%
% OUTPUT:   Tx_RGB: horizontal shifts, a 40x3 matrix
%           Ty_RGB: vertical shifts, a 40x3 matrix
%
% NOTE: _Tx_RGB(1,:) = Ty_RGB(1,:) = (0 0 0) by definition.
%       _Tx_RGB(20,2) is the horizontal shift of the Green layer of the
%       20th image relatively to the Green layer of the firs image.
%
%
% OUTLINE OF THE ALGORITHM:
%
% 1.The first step is to compute the continuous moments m_00, m_01 and m_10
% of each low-resolution image using the .mat file called:
% PolynomialReproduction_coef.mat. This file contains three matrices
% 'Coef_0_0', 'Coef_1_0' and 'Coef_0_1' used to calculate the continuous
% moments.
%
% 2.The second step consists in calculating the barycenters of the Red,
% Green and Blue layers of the low-resolution images.
%
% 3.By computing the difference between the barycenters of corresponding
% layers between two images, the horizontal and vertical shifts can be
% retrieved for each layer.
%
%
% Author:   Loic Baboulaz
% Date:     August 2006
%
% Imperial College London
% *************************************************************************
% Load the coefficients for polynomial reproduction

load('PolynomialReproduction_coef.mat', 'Coef_0_0', 'Coef_1_0', 'Coef_0_1');
% -------- include your code here -----------
% parameters
N = 40; % number of cameras
n = 3; % number of layers
th = 0.28; % noise threshold 
x = zeros(N, n);
y = zeros(N, n);
Tx_RGB = zeros(N, n);
Ty_RGB = zeros(N, n);
% load low-resolution images
for i = 1: N    
    data = double(imread(sprintf('LR_Tiger_%.2d.tif', i)))/ 255; % obtain and rescale samples   
    for j = 1: n
        % reduce noise by comparing the value with threshold 
        store = data(:, :, j);
        if j == 1 % red layer
            store(store < th) = 0;
            data(:, :, j) = store;
        elseif j == 2 % green layer
            store(store < th) = 0;
            data(:, :, j) = store;
        else % blue layer
            store(store < th) = 0;
            data(:, :, j) = store;
        end
        m_0_0 = sum(sum(Coef_0_0 .* data(:, :, j)));  % moments
        m_0_1 = sum(sum(Coef_0_1 .* data(:, :, j)));
        m_1_0 = sum(sum(Coef_1_0 .* data(:, :, j)));       
        x(i, j) = m_1_0 / m_0_0; % centralize
        y(i, j) = m_0_1 / m_0_0;
    end
    % the first figure is the reference and calculate the shift 
    Tx_RGB(i, :) = x(i, :) - x(1, :);
    Ty_RGB(i, :) = y(i, :) - y(1, :);
end
end