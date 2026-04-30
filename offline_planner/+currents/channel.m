function [v_c] = channel(~, ~, y, y_center, channel_width, intensity)
% Compute channel current in the x-direction.
%
% Inputs:
%   t             - time [s] (unused, time-invariant)
%   x             - x positions [m] (unused)
%   y             - y positions [m]
%   y_center      - center of the channel [m]
%   channel_width - total width of the channel [m]
%   intensity     - peak flow speed [m/s]
%
% Outputs:
%   v_c - current velocity (2 x n) [m/s]

y_normalized = abs(y - y_center) / (channel_width / 2);
v_x          = intensity * max(0, 1 - y_normalized.^2);
v_c          = [v_x; zeros(size(v_x))];
end
