function [z] = interpolate(y, tau, tauf)
% Interpolate collocation data using Lagrange polynomials.
%
% Inputs:
%   y    - values at collocation points (n x N+1)
%   tau  - query times [s]
%   tauf - final time [s]
%
% Outputs:
%   z - interpolated values (n x numel(tau))

% extract sizes
[n, M] = size(y);
N      = M - 1;

% compute legendre polynomial and its derivative
[L, Ldot] = lpm.compute_legendre_polynomial(N);

% compute collocation points
col = lpm.compute_collocation_points(N);

% form scaled time array
t = 2 * tau / tauf - 1;

% interpolate
z = zeros(n, numel(t));

for kk=1:numel(col)
    % compute value of lagrange polynomial
    scl = 1 / (N * (N + 1)) / polyval(L, col(kk));
    phi = scl * (t.^2 - 1) .* polyval(Ldot, t) ./ (t - col(kk));

    % account for collocation points (division by zero above)
    cidx      = isnan(phi) | isinf(phi);
    phi(cidx) = 1;

    % add to interpolation
    z = z + y(:,kk) * phi;
end
end