function [H_anom, R_anom] = anomaly_measurement(pE, pN, g_field, m_field, origin)
% Compute the anomaly measurement Jacobian and noise covariance at vehicle position.
%
% Noise standard deviations from:
%   Brodovsky & Dames, "Navigation in GNSS-Denied Environments Using
%   MEMS-Grade Sensors and Geophysical Anomalies: A UKF Approach"
%   sigma_g = 48.93 mGal, sigma_m = 41,303 nT
%
% Inputs:
%   pE      - east position [m]
%   pN      - north position [m]
%   g_field - function handle mapping (x, y) to (1 x n) gravity anomaly [mGal]
%   m_field - function handle mapping (x, y) to (1 x n) magnetic anomaly [nT]
%   origin  - geographic origin [lat0, lon0] [deg]
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

sigma_g = 48.93;
sigma_m = 41303;
R_anom  = diag([sigma_g, sigma_m] .^ 2);
end
