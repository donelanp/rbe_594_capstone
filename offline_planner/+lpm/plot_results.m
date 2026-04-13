function [] = plot_results(a, b, tauf, x0, xf, dynfun)
% Plot optimized state and control trajectories.
%
% Inputs:
%   a      - state at collocation points (n x N+1)
%   b      - control at collocation points (m x N+1)
%   tauf   - final time [s]
%   x0     - initial state (n x 1)
%   xf     - final state (n x 1)
%   dynfun - dynamics function handle f(tau, x, u)

% dimenions of system state and control input
n = size(a, 1);
m = size(b, 1);

% interpolate numerically computed state and control inputs
tau_i = linspace(0, tauf);
x_i   = lpm.interpolate(a, tau_i, tauf);
u_i   = lpm.interpolate(b, tau_i, tauf);

% simulate the system dynamics with the optimized control inputs
[tau_s, x_s] = lpm.simulate(x0, b, tauf, dynfun);

% plot system state
figure('Theme', 'light');
axn  = gobjects(n, 1);
tile = tiledlayout(n, 1, 'TileIndexing', 'columnmajor');

for ii=1:n
    axn(ii) = nexttile(tile);
    hold on;
    plot(tau_i, x_i(ii,:), 'LineWidth', 2, 'DisplayName', 'interpolated');
    plot(tau_s, x_s(:,ii), '--', 'LineWidth', 2, 'DisplayName', 'simulated');
    plot([0 tauf], [x0(ii) xf(ii)], 'k^', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'boundary conditions');
    hold off;
    ylabel(sprintf('x_%d', ii));

    if ii == 1
        title('state vs time');
        legend();
    elseif ii == n
        xlabel('time [s]');
    end
end

linkaxes(axn, 'x');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold');

% plot control input
figure('Theme', 'light');
axm  = gobjects(m, 1);
tile = tiledlayout(m, 1, 'TileIndexing', 'columnmajor');

for ii=1:m
    axm(ii) = nexttile(tile);
    plot(tau_i, u_i(ii,:), 'LineWidth', 2);
    ylabel(sprintf('u_%d', ii));

    if ii == 1
        title('control input vs time');
    elseif ii == m
        xlabel('time [s]');
    end
end

linkaxes(axm, 'x');
set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold');
end