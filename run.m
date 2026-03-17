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

% initial covariance P0 and its lower cholesky factor L0 (P0 = L0*L0')
sigma_p   = 10;                    % initial position std [m]
sigma_v   = 0.1;                   % initial velocity std [m/s]
sigma_psi = deg2rad(0.01);         % initial heading std [rad]
sigma_ba  = 30e-6 * g;             % initial accel bias std [m/s^2]
sigma_bg  = deg2rad(0.001) / 3600; % initial gyro bias std [rad/s]
P0        = diag([sigma_p, sigma_p, sigma_v, sigma_v, sigma_psi, sigma_ba, sigma_ba, sigma_bg] .^ 2);
L0        = chol(P0, 'lower');
L0_vec    = L_to_lvec(L0);

% non-dimensionalization of vehicle states
L_char = norm(xf - x0(1:2)); % characteristic length [m]
T_char = L_char / max_speed; % characteristic time [s]

x0_nd            = [x0(1)/L_char; x0(2)/L_char; x0(3)];
xf_nd            = [xf(1)/L_char; xf(2)/L_char];
current_field_nd = @(t_nd, x_nd, y_nd) current_field(T_char * t_nd, L_char * x_nd, L_char * y_nd) / max_speed;
wheelbase_nd     = wheelbase / L_char;
max_speed_nd     = 1;
max_steer_nd     = max_steer;
max_accel_nd     = max_accel * T_char / max_speed;
max_srate_nd     = max_srate * T_char;

% non-dimensionalization of covariance
D_s_vec   = [L_char, L_char, max_speed, max_speed, 1, max_accel, max_accel, max_srate];
D_s_inv   = diag(1 ./ D_s_vec);
L0_nd     = D_s_inv * L0;
L0_vec_nd = L_to_lvec(L0_nd);
Q_nd      = T_char * D_s_inv * Q * D_s_inv;

% cost function
cost = @(z) compute_cost(z, N);

% constraint function (dynamics and boundary conditions)
% full state: x = [pE_nd; pN_nd; theta; L_nd_vec (36 elements)] = 39 states
% z = [x (39*(N+1)); u (2*(N+1)); tauf (1)]
constraints = @(z) nonlcon(z, x0_nd, xf_nd, L0_vec_nd, current_field_nd, wheelbase_nd, max_speed_nd, max_steer_nd, max_accel_nd, max_srate_nd, Q_nd, N);

% initial guess (straight line, constant controls, constant covariance)
z0 = initialize_guess(N, x0_nd, xf_nd, L0_vec_nd);

% decision variable bounds
diag_ind = [1, 9, 16, 22, 27, 31, 34, 36];
lb_state = [-Inf; -Inf; -pi; -Inf(36, 1)];
ub_state = [ Inf;  Inf;  pi;  Inf(36, 1)];
for idx = diag_ind
    lb_state(3 + idx) = eps();
end
lb = [repmat(lb_state, N+1, 1); repmat([-1; -1], N+1, 1); 0  ];
ub = [repmat(ub_state, N+1, 1); repmat([ 1;  1], N+1, 1); Inf];

% optimize
zopt = fmincon(cost, z0, [], [], [], [], lb, ub, constraints);

% simulate vehicle based on optimal control inputs
[x_nd, u, tauf_nd] = extract_state(zopt, N);
tauf               = tauf_nd * T_char;
x_uuv              = x_nd(1:3, :) .* [L_char; L_char; 1];
x_cov              = x_nd(4:end,:);

dynfun_uuv_dim = @(tau, x, u) uuv_dynamics(tau, x, u, current_field, wheelbase, max_speed, max_steer);
[t_s, x_s]     = lpm.simulate(x0, u, tauf, dynfun_uuv_dim);
t_i            = linspace(0, tauf, 100);
x_i            = lpm.interpolate(x_uuv, t_i, tauf);

% plot state and control input vs time
lpm.plot_results(x_uuv, u, tauf, x0, [xf; nan], dynfun_uuv_dim);

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

    % current speed
    C     = current_field(repmat(t_i(k), 1, numel(X1)), X1(:)', X2(:)');
    speed = reshape(sqrt(C(1,:).^2 + C(2,:).^2), size(X1));
    pcolor(X1 / 1000, X2 / 1000, speed, 'EdgeColor', 'none', 'FaceAlpha', 0.7, 'FaceColor', 'interp', 'HandleVisibility', 'off');

    % current direction
    Cq = current_field(repmat(t_i(k), 1, numel(X1)), X1(:)', X2(:)');
    quiver(X1(:)' / 1000, X2(:)' / 1000, Cq(1,:), Cq(2,:), 'k', 'LineWidth', 0.5, 'HandleVisibility', 'off');

    % vehicle position
    idx_s = t_s <= t_i(k);
    plot(x0(1) / 1000, x0(2) / 1000, 'k^', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'initial position');
    plot(xf(1) / 1000, xf(2) / 1000, 'kd', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'goal position');
    plot(x_i(1,1:k) / 1000, x_i(2,1:k) / 1000, 'g-', 'LineWidth', 2, 'DisplayName', 'optimal trajectory');
    plot(x_s(idx_s,1) / 1000, x_s(idx_s,2) / 1000, 'r--', 'LineWidth', 2, 'DisplayName', 'simulated trajectory');

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

%% compute cost function to minimize control input
function [J] = compute_cost(z, N)
[~, u, tauf] = extract_state(z, N);

w = lpm.compute_quadrature_weights(N);
L = 0.5 * sum(u.^2, 1);
J = (tauf / 2) * sum(w .* L);
end

%% initialize decision variables
function [z0] = initialize_guess(N, x0, xf, L0_vec_nd)
% non-dimensional final time guess (~2x straight-line travel time)
tauf = 2;

% heading toward goal (atan2 of east/north for compass convention)
theta0 = atan2(xf(1) - x0(1), xf(2) - x0(2));

% straight line interpolation for position, constant heading, constant covariance
tau = lpm.compute_collocation_points(N);
s   = 0.5 * (tau + 1);
xy  = x0(1:2) * (1 - s) + xf * s;
x   = [xy; theta0 * ones(1, N+1); repmat(L0_vec_nd, 1, N+1)];

% constant half speed, zero steering (normalized)
u = [0.5 * ones(1, N+1); zeros(1, N+1)];

z0 = [x(:); u(:); tauf];
end

%% extract state variables from decision vector
function [x, u, tauf] = extract_state(z, N)
ix = 39 * (N + 1);
iu = ix + 2 * (N + 1);

x    = reshape(z(1:ix), 39, []);
u    = reshape(z(ix+1:iu), 2, []);
tauf = z(end);
end

%% compute constraints
function [c, ceq] = nonlcon(z, x0, xf, L0_vec, current_field, wheelbase, max_speed, max_steer, max_accel, max_srate, Q, N)
% extract state variables
[x, u, tauf] = extract_state(z, N);

% compute collocation points
tau = 0.5 * tauf * (lpm.compute_collocation_points(N) + 1);

% compute differentiation matrix
D = lpm.compute_differentiation_matrix(N);

% compute inequality constraints due to acceleration and steering rate
% limits
u_dot      = (2 / tauf) * u * D';
accel      = u_dot(1,:) * max_speed;
steer_rate = u_dot(2,:) * max_steer;

c = [ accel - max_accel;
     -accel - max_accel;
      steer_rate - max_srate;
     -steer_rate - max_srate]';
c = c(:);

% compute equality constraints due to system dynamics
ceq       = zeros(39 * (N + 1) + 41, 1);
ic        = 39 * (N + 1);
ceq(1:ic) = reshape(x * D' - 0.5 * tauf * system_dynamics(tau, x, u, current_field, wheelbase, max_speed, max_steer, max_accel, max_srate, Q), [], 1);

% compute equality constraints due to boundary conditions
ceq(ic+1:ic+3)   = x(1:3,   1) - x0;     % initial vehicle state
ceq(ic+4:ic+39)  = x(4:end, 1) - L0_vec; % initial error covariance
ceq(ic+40:ic+41) = x(1:2,  end) - xf;    % final vehicle state
end

%% uuv dynamics of kinematic vehicle with current
function [xdot] = uuv_dynamics(tau, x, u, current_field, wheelbase, max_speed, max_steer)
theta = x(3,:);
v     = u(1,:) * max_speed;
delta = u(2,:) * max_steer;
c     = current_field(tau, x(1,:), x(2,:));

xdot = [v .* sin(theta) + c(1,:);
        v .* cos(theta) + c(2,:);
        (v / wheelbase) .* tan(delta)];
end

%% covariance dynamics
function [Ldot_vec] = cov_dynamics(psi, x_cov, max_accel_nd, max_srate_nd, Q)
n_pts    = size(x_cov, 2);
Ldot_vec = zeros(36, n_pts);

for k = 1:n_pts
    % form covariance matrix from cholesky factorization
    L = lvec_to_L(x_cov(:, k));
    P = L * L';

    % form non-dimensional linearized dynamics Jacobian
    cp = cos(psi(k));
    sp = sin(psi(k));
    F  = zeros(8);

    F(1, 3) = 1;
    F(2, 4) = 1;
    F(3, 6) = max_accel_nd * cp;
    F(3, 7) = -max_accel_nd * sp;
    F(4, 6) = max_accel_nd * sp;
    F(4, 7) = max_accel_nd * cp;
    F(5, 8) = max_srate_nd;

    % compute covariance derivative using differential lyapunov function
    Pdot = F * P + P * F' + Q;

    % compute derivative of cholesky factorization
    M              = L \ Pdot / L';
    S              = tril(M, -1) + 0.5 * diag(diag(M));
    Ldot           = L * S;
    Ldot_vec(:, k) = L_to_lvec(Ldot);
end
end

%% full system dynamics
function [xdot] = system_dynamics(tau, x, u, current_field, wheelbase, max_speed, max_steer, max_accel_nd, max_srate_nd, Q)
xdot = [uuv_dynamics(tau, x(1:3,:), u, current_field, wheelbase, max_speed, max_steer);
        cov_dynamics(x(3,:), x(4:end,:), max_accel_nd, max_srate_nd, Q)];
end

%% lower triangular matrix <==> lower triangular element vector
function [l_vec] = L_to_lvec(L)
l_vec = L(tril(true(8,8)));
end

function [L] = lvec_to_L(l_vec)
L                  = zeros(8,8);
L(tril(true(8,8))) = l_vec;
end
