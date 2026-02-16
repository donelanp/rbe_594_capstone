function [tau, x] = simulate(x0, b, tauf, dynfun)
% Simulate system dynamics with interpolated control inputs.
%
% Inputs:
%   x0     - initial state (n x 1)
%   b      - control inputs at collocation points (m x N+1)
%   tauf   - final time [s]
%   dynfun - dynamics function handle f(tau, x, u)
%
% Outputs:
%   tau - time vector [s]
%   x   - state trajectory (numel(tau) x n)

% wrap the dynamic function to use interpolated control input
odefun = @(tau, x) dynfun(tau, x, lpm.interpolate(b, tau, tauf));

% simulate the system
[tau, x] = ode45(odefun, [0 tauf], x0);
end