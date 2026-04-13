function [v_c] = oscillating_uniform(t, ~, ~, intensity, frequency)
% Compute spatially uniform oscillating current in the x-direction.
%
% Inputs:
%   t         - time [s]
%   x         - x positions [m] (unused, spatially uniform)
%   y         - y positions [m] (unused, spatially uniform)
%   intensity - peak flow speed [m/s]
%   frequency - oscillation frequency [Hz]
%
% Outputs:
%   v_c - current velocity (2 x n) [m/s]

n   = numel(t);
v_x = intensity * cos(2 * pi * frequency * t);
v_c = [v_x; zeros(1, n)];
end
