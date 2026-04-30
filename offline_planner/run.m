% number of collocation intervals
N = 23;

% boundary conditions [m, m, rad]
% orientation = 0 is north, orientation = pi/2 is east
% final orientation is free
x0 = [-500000; -200000; 0];
xf = [-350000; -120000];

% current model selection
current_model = 'hycom';
switch current_model
    case 'zero'
        current_field = @(t, x, y) zeros(2, numel(x));
    case 'channel'
        y_center      = -85000;
        channel_width = 230000;
        intensity     = 1.0;
        current_field = @(t, x, y) currents.channel(t, x, y, y_center, channel_width, intensity);
    case 'shear'
        intensity     = 2.5e-5;
        current_field = @(t, x, y) currents.shear(t, x, y, intensity);
    case 'vortex'
        x_center      = -225000;
        y_center      = -85000;
        intensity     = 40000;
        current_field = @(t, x, y) currents.vortex(t, x, y, x_center, y_center, intensity);
    case 'oscillating_uniform'
        intensity     = 1.0;
        frequency     = 1.5e-5;
        current_field = @(t, x, y) currents.oscillating_uniform(t, x, y, intensity, frequency);
    case 'hycom'
        origin        = [-37.5, 22.5];
        region        = [-660000, 660000, -500000, 500000];
        current_field = @(t, x, y) currents.hycom(t, x, y, origin, region);
    otherwise
        error('Unknown current model: %s', current_model);
end

% dimensional vehicle parameters
wheelbase = 1;      % wheelbase [m] (bicycle kinematic model)
max_speed = 2;      % maximum vehicle speed applied by motor [m/s]
max_steer = pi / 3; % maximum steering angle [rad]
max_accel = 1;      % maximum acceleration [m/s^2]
max_srate = pi / 4; % maximum steering rate [rad/s]

% dimensional inertial navigation process noise covariance matrix
% for error state [dp_N, dp_E, dv_N, dv_E, dpsi, bax, bay, bgz] assuming 
% aviation grade IMU
g     = 9.80665;
S_ra  = (40e-6 * g) ^ 2;                        % accelerometer random walk [m^2/s^3]
S_rg  = (deg2rad(0.24) / 3600) ^ 2;             % gyroscope random walk [rad^2/s]
S_bad = (5.585216e-6 * g) ^ 2;                  % accelerometer bias instability [m^2/s^5]
S_bgd = (deg2rad(10.501290) / 3600 / 3600) ^ 2; % gyroscope bias instability [rad^2/s^3]
Q     = diag([0, 0, S_ra, S_ra, S_rg, S_bad, S_bad, S_bgd]);

% DVL and compass measurement model for error state
% [dp_N, dp_E, dv_N, dv_E, dpsi, bax, bay, bgz]
sigma_dvl     = 0.01;         % DVL velocity noise std [m/s]
sigma_compass = deg2rad(0.5); % compass heading noise std [rad]

H = [zeros(2, 2), -eye(2, 2), zeros(2, 4); zeros(1, 4), -1, zeros(1,3)];
R = diag([sigma_dvl, sigma_dvl, sigma_compass] .^ 2);

% initial covariance P0 and its lower cholesky factor L0 (P0 = L0*L0')
sigma_p   = 10;                    % initial position std [m]
sigma_v   = 0.1;                   % initial velocity std [m/s]
sigma_psi = deg2rad(0.01);         % initial heading std [rad]
sigma_ba  = 30e-6 * g;             % initial accel bias std [m/s^2]
sigma_bg  = deg2rad(0.001) / 3600; % initial gyro bias std [rad/s]
P0        = diag([sigma_p, sigma_p, sigma_v, sigma_v, sigma_psi, sigma_ba, sigma_ba, sigma_bg] .^ 2);
Lp0       = chol(P0(1:2, 1:2), 'lower');
Lp0_vec   = utils.L_to_lvec(Lp0);

% steady-state velocity error covariance
Lvss_vec = nav.vel_cov_steady_state(Q, H, R);

% parameters used for non-dimensionalization
L_char = norm(xf - x0(1:2)); % characteristic length [m]
T_char = L_char / max_speed; % characteristic time [s]
Gamma  = [L_char; L_char; pi; L_char; L_char; L_char];

% cost function weight (alpha = 0: control effort only, alpha = 1: covariance only)
alpha = 1.0;
cost  = @(z) compute_cost(z, N, alpha);

% constraint function (dynamics and boundary conditions)
% x_nd = [pE_nd (1); pN_nd (1); theta_nd (1); Lp_vec_nd (3)]
% u_nd = [v_nd (1); δ_nd (1)]
% z    = [x_nd (6*(N+1)); u_nd (2*(N+1)); tauf_nd (1)]
constraints = @(z) nonlcon(z, x0, xf, Lp0_vec, current_field, wheelbase, max_speed, max_steer, max_accel, max_srate, T_char, Gamma, Lvss_vec, N);

% initial guess
z0 = initialize_guess(N, x0, xf, Lp0_vec, T_char, Lvss_vec, Gamma);

% decision variable bounds
lb_state = [-Inf; -Inf; -1; eps(); -Inf; eps()];
ub_state = [ Inf;  Inf;  1;  Inf;  Inf;  Inf];
lb = [repmat(lb_state, N+1, 1); repmat([0; -1], N+1, 1); 0  ];
ub = [repmat(ub_state, N+1, 1); repmat([ 1;  1], N+1, 1); Inf];

% optimize
opts = optimoptions('fmincon', ...
    'Algorithm',              'sqp', ...
    'Display',                'iter-detailed', ...
    'MaxIterations',          10000, ...
    'MaxFunctionEvaluations', 10000);
zopt = fmincon(cost, z0, [], [], [], [], lb, ub, constraints, opts);

% simulate vehicle based on optimal control inputs
[x_nd, u_nd, tauf_nd] = extract_state(zopt, N);
tauf                  = tauf_nd * T_char;
u                     = u_nd .* [max_speed; max_steer];
x_uuv                 = x_nd(1:3, :) .* [L_char; L_char; pi];

dynfun_uuv = @(tau, x, u) nav.uuv_dynamics(tau, x, u, current_field, wheelbase);
[t_s, x_s] = lpm.simulate(x0, u, tauf, dynfun_uuv);
t_i        = linspace(0, tauf, 100);
x_i        = lpm.interpolate(x_uuv, t_i, tauf);

% plot state and control input vs time
lpm.plot_results(x_uuv, u, tauf, x0, [xf; nan], dynfun_uuv);

% sample spatial region to plot current magnitude and direction
margin           = 0.2;
[x1_min, x1_max] = bounds(x_uuv(1,:));
[x2_min, x2_max] = bounds(x_uuv(2,:));
x1_pad           = margin * (x1_max - x1_min);
x2_pad           = margin * (x2_max - x2_min);
x1               = linspace(x1_min - x1_pad, x1_max + x1_pad, 41);
x2               = linspace(x2_min - x2_pad, x2_max + x2_pad, 41);
[X1, X2]         = meshgrid(x1, x2);

% compute max current speed across all frames for fixed color scale
max_cspeed = eps();
for k = 1:numel(t_i)
    C          = current_field(repmat(t_i(k), 1, numel(X1)), X1(:)', X2(:)');
    max_cspeed = max(max_cspeed, max(vecnorm(C, 2, 1)));
end

% plot vehicle trajectory over time
fig           = figure('Theme', 'light', 'Color', 'w', 'Position', [100 100 1024 512]);
vid           = VideoWriter(['trajectory_' current_model '.mp4'], 'MPEG-4');
vid.FrameRate = 30;
open(vid);
for k = 1:numel(t_i)
    clf;
    hold on;

    % current speed and direction
    C     = current_field(repmat(t_i(k), 1, numel(X1)), X1(:)', X2(:)');
    speed = reshape(vecnorm(C, 2, 1), size(X1));
    pcolor(X1 / 1000, X2 / 1000, speed, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'FaceColor', 'interp', 'HandleVisibility', 'off');
    quiver(X1(:)' / 1000, X2(:)' / 1000, C(1,:), C(2,:), 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % vehicle position
    idx_s = t_s <= t_i(k);
    plot(x0(1) / 1000, x0(2) / 1000, 'k^', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'initial position');
    plot(xf(1) / 1000, xf(2) / 1000, 'kd', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'goal position');
    plot(x_i(1,1:k) / 1000, x_i(2,1:k) / 1000, 'g-', 'LineWidth', 2, 'DisplayName', 'optimal trajectory');
    plot(x_s(idx_s,1) / 1000, x_s(idx_s,2) / 1000, 'r--', 'LineWidth', 2, 'DisplayName', 'reconstructed trajectory');

    % housekeeping
    hold off;

    cb              = colorbar;
    cb.Label.String = 'Current Speed [m/s]';
    xlabel('East [km]');
    ylabel('North [km]');
    title(sprintf('vehicle trajectory through %s current (t = %.1f hrs)', current_model, t_i(k) / 3600), 'Interpreter', 'none');
    legend('Location', 'northwest');
    axis equal;
    xlim([(x1_min - x1_pad) / 1000, (x1_max + x1_pad) / 1000]);
    ylim([(x2_min - x2_pad) / 1000, (x2_max + x2_pad) / 1000]);
    clim([0, max_cspeed]);

    set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
    set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold');

    drawnow;
    writeVideo(vid, getframe(gcf));
end
close(vid);

% plot total position covariance vs time
Lp_i    = lpm.interpolate(x_nd(4:6,:) .* Gamma(4:6), t_i, tauf);
trace_i = Lp_i(1,:).^2 + Lp_i(2,:).^2 + Lp_i(3,:).^2;

theta_s        = @(t) interp1(t_s, x_s(:,3), t, 'linear', 'extrap');
L0_vec         = utils.L_to_lvec(chol(P0, 'lower'));
[t_rec, L_rec] = ode15s(@(t, L_vec) nav.cov_dynamics(-theta_s(t), L_vec, Q, H, R), [0, tauf], L0_vec);
trace_rec      = zeros(size(t_rec));

for k = 1:numel(t_rec)
    Lk           = utils.lvec_to_L(L_rec(k,:)');
    Pk           = Lk * Lk';
    trace_rec(k) = Pk(1,1) + Pk(2,2);
end

figure('Theme', 'light', 'Color', 'w');
hold on;
plot(t_i / 3600, trace_i / 1e6, 'g-', 'LineWidth', 2, 'DisplayName', 'optimal covariance');
plot(t_rec / 3600, trace_rec / 1e6, 'r--', 'LineWidth', 2, 'DisplayName', 'reconstructed covariance');
hold off;
xlabel('Time [hr]');
ylabel('\sigma_{pos}^2 [km^2]');
title('position covariance vs time');
legend('Location', 'northwest');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold');

%% compute cost function
function [J] = compute_cost(z, N, alpha)
[x_nd, u_nd, tauf_nd] = extract_state(z, N);

w         = lpm.compute_quadrature_weights(N);
ctrl_cost = 0.5 * sum(u_nd.^2, 1);
cov_cost  = x_nd(4,:).^2 + x_nd(5,:).^2 + x_nd(6,:).^2;
L         = (1 - alpha) * ctrl_cost + alpha * cov_cost;
J         = (tauf_nd / 2) * sum(w .* L);
end

%% initialize decision variables
function [z0] = initialize_guess(N, x0, xf, Lp0_vec, T_char, Lvss_vec, Gamma)
% straight-line travel at full-tilt
u_nd = [ones(1, N+1); zeros(1, N+1)];

% duration for straight-line travel at full-tilt
tauf    = T_char;
tauf_nd = tauf / T_char;

% straight line interpolation for position along constant heading
tau = lpm.compute_collocation_points(N);
s   = 0.5 * (tau + 1);
xy  = x0(1:2) * (1 - s) + xf * s;

% heading toward goal (atan2 of east/north for compass convention)
theta0 = atan2(xf(1) - x0(1), xf(2) - x0(2)) * ones(1, N+1);

% propagate covariance along straight-line constant heading trajectory
t          = 0.5 * tauf * (tau + 1);
[~, Lp_vec] = ode45(@(t, Lp_vec) nav.pos_cov_dynamics(Lp_vec, Lvss_vec), t, Lp0_vec);

x_nd = [xy; theta0; Lp_vec'] ./ Gamma;
z0   = [x_nd(:); u_nd(:); tauf_nd];
end

%% extract state variables from decision vector
function [x_nd, u_nd, tauf_nd] = extract_state(z, N)
ix = 6 * (N + 1);
iu = ix + 2 * (N + 1);

x_nd    = reshape(z(1:ix), 6, []);
u_nd    = reshape(z(ix+1:iu), 2, []);
tauf_nd = z(end);
end

%% compute constraints
function [c, ceq] = nonlcon(z, x0, xf, Lp0_vec, current_field, wheelbase, max_speed, max_steer, max_accel, max_srate, T_char, Gamma, Lvss_vec, N)
% extract state variables
[x_nd, u_nd, tauf_nd] = extract_state(z, N);

% re-dimensionalize
u    = u_nd .* [max_speed; max_steer];
x    = x_nd .* Gamma;
tauf = tauf_nd * T_char;

% compute collocation points
tau = 0.5 * tauf * (lpm.compute_collocation_points(N) + 1);

% compute differentiation matrix
D = lpm.compute_differentiation_matrix(N);

% compute inequality constraints due to acceleration and steering rate
% limits
u_dot      = (2 / tauf) * u * D';
accel      = u_dot(1,:);
steer_rate = u_dot(2,:);

c = [ (accel      - max_accel) / max_accel;
     -(accel      + max_accel) / max_accel;
      (steer_rate - max_srate) / max_srate;
     -(steer_rate + max_srate) / max_srate]';
c = c(:);

% compute equality constraints due to system dynamics
ceq       = zeros(6 * (N + 1) + 8, 1);
ic        = 6 * (N + 1);
ceq(1:ic) = reshape((x * D' - 0.5 * tauf * nav.system_dynamics(tau, x, u, current_field, wheelbase, Lvss_vec)) ./ Gamma, [], 1);

% compute equality constraints due to boundary conditions
ceq(ic+1:ic+3) = (x(1:3, 1) - x0) ./ Gamma(1:3);          % initial vehicle state
ceq(ic+4:ic+6) = (x(4:end, 1) - Lp0_vec) ./ Gamma(4:end); % initial position covariance
ceq(ic+7:ic+8) = (x(1:2, end) - xf) ./ Gamma(1:2);        % final vehicle state
end
