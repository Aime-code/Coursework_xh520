function [phi_T] = bspline(period, p)
% obtain the sample of bspline of order N 
% Input:
% period: sampling period
% p: the order of bspline
% Output:
% phi_T: sample of bspline of a given order

boxFun = ones(1, period); % b-spline of order 0 is the box function
phiT = boxFun; 
if p == 0
    phi_T = boxFun;
else
    for i = 1: p
        if mod(i, 2)
            phi_T = [0 conv(boxFun, phiT)] / period; % centralize
        else
            phi_T = [conv(boxFun, phiT) 0] / period;     
        end
        phiT = phi_T;
    end
end
end
