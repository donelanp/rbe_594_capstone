function [v_c] = shear(x, y, intensity)
% Compute linear shear current in the x-direction.
%
% Inputs:
%   x         - x positions [m]
%   y         - y positions [m]
%   intensity - shear rate [1/s]
%
% Outputs:
%   v_c - current velocity (2 x n) [m/s]

v_c = [intensity * y; zeros(size(y))];
end
