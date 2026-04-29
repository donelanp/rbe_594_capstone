function [H_anom, R_anom] = anomaly_measurement(pE, pN, g_field, m_field, origin, region)
% Compute the anomaly measurement Jacobian and noise covariance at vehicle position.
%
% Noise standard deviations approximated as the standard deviation of the
% map's high-frequency residual (field - Gaussian-smoothed field) sampled
% over the scenario region. The idea is to capture map interpolation error
% and sensor noise into a single value.
%
% Inputs:
%   pE      - east position [m]
%   pN      - north position [m]
%   g_field - function handle mapping (x, y) to (1 x n) gravity anomaly [mGal]
%   m_field - function handle mapping (x, y) to (1 x n) magnetic anomaly [nT]
%   origin  - geographic origin [lat0, lon0] [deg]
%   region  - scenario bounds [east_min, east_max, north_min, north_max] [m]
%
% Outputs:
%   H_anom - measurement Jacobian (2 x 2)
%   R_anom - measurement noise covariance matrix (2 x 2)

% data grid spacing at origin latitude [m]
m_per_arcmin = pi / 180 / 60 * 6371000 * cos(deg2rad(origin(1)));
dp_g         = 1 * m_per_arcmin;
dp_m         = 3 * m_per_arcmin;

% gradient approximations using central difference
H_anom = -[(g_field(pE, pN+dp_g) - g_field(pE, pN-dp_g)) / (2*dp_g), ...
          (g_field(pE+dp_g, pN) - g_field(pE-dp_g, pN)) / (2*dp_g);
          (m_field(pE, pN+dp_m) - m_field(pE, pN-dp_m)) / (2*dp_m), ...
          (m_field(pE+dp_m, pN) - m_field(pE-dp_m, pN)) / (2*dp_m)];

persistent sigma_g sigma_m
if isempty(sigma_g)
    [sigma_g, sigma_m] = estimate_sigma(g_field, m_field, region);
end

R_anom = diag([sigma_g, sigma_m] .^ 2);
end

function [sigma_g, sigma_m] = estimate_sigma(g_field, m_field, region)
nx     = 201;
ny     = 201;
x      = linspace(region(1), region(2), nx);
y      = linspace(region(3), region(4), ny);
[X, Y] = meshgrid(x, y);

G        = reshape(g_field(X(:)', Y(:)'), ny, nx);
M        = reshape(m_field(X(:)', Y(:)'), ny, nx);
G_smooth = imgaussfilt(G, 2);
M_smooth = imgaussfilt(M, 2);

sigma_g = std(G(:) - G_smooth(:));
sigma_m = std(M(:) - M_smooth(:));
end
