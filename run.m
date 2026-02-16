% number of collocation intervals
N = 15;

% boundary conditions [m, m, rad]
% orientation = 0 is north, orientation = pi/2 is east
% final orientation is free
x0 = [0; 0; 0];
xf = [100; 50];

% current model selection
current_model = 'channel';
switch current_model
    case 'channel'
        y_center      = 25;
        channel_width = 50;
        intensity     = 1.0;
        current_field = @(x, y) currents.channel(x, y, y_center, channel_width, intensity);
    case 'shear'
        intensity     = 0.1;
        current_field = @(x, y) currents.shear(x, y, intensity);
    case 'vortex'
        x_center      = 50;
        y_center      = 25;
        intensity     = 5;
        current_field = @(x, y) currents.vortex(x, y, x_center, y_center, intensity);
    otherwise
        error('Unknown current model: %s', current_model);
end

% vehicle parameters
wheelbase = 1;      % wheelbase [m] (bicycle kinematic model)
max_speed = 2;      % maximum vehicle speed applied by motor [m/s]
max_steer = pi / 3; % maximum steering angle [rad]

% cost function
% z = [x (3*(N+1)); u (2*(N+1)); tauf (1)]
cost = @(z) compute_cost(z, N);

% constraint function (dynamics and boundary conditions)
% z = [x (3*(N+1)); u (2*(N+1)); tauf (1)]
constraints = @(z) nonlcon(z, x0, xf, current_field, wheelbase, N);

% initial guess (straight line, constant controls)
z0 = initialize_guess(N, x0, xf);

% decision variable bounds
lb = [repmat([-Inf; -Inf; -pi], N+1, 1); repmat([-max_speed; -max_steer], N+1, 1); 0];
ub = [repmat([ Inf;  Inf;  pi], N+1, 1); repmat([ max_speed;  max_steer], N+1, 1); Inf];

% optimize
zopt = fmincon(cost, z0, [], [], [], [], lb, ub, constraints);

% plot results
[x, u, tauf] = extract_state(zopt, N);
lpm.plot_results(x, u, tauf, x0, [xf; nan], @(tau, x, u) system_dynamics(tau, x, u, current_field, wheelbase));

figure('Theme', 'light');
hold on;

margin               = 0.2;
[x1_min, x1_max]     = bounds(x(1,:));
[x2_min, x2_max]     = bounds(x(2,:));
x1_pad               = margin * (x1_max - x1_min);
x2_pad               = margin * (x2_max - x2_min);

x1                   = linspace(x1_min - x1_pad, x1_max + x1_pad, 25);
x2                   = linspace(x2_min - x2_pad, x2_max + x2_pad, 25);
[X1, X2]             = meshgrid(x1, x2);
X1                   = reshape(X1, 1, []);
X2                   = reshape(X2, 1, []);
C                    = current_field(X1, X2);
C1                   = C(1,:);
C2                   = C(2,:);
quiver(X1, X2, C1, C2, 'DisplayName', 'current');

[~, x_s] = lpm.simulate(x0, u, tauf, @(tau, x, u) system_dynamics(tau, x, u, current_field, wheelbase));

plot(x0(1), x0(2), 'k^', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'initial position');
plot(xf(1), xf(2), 'kd', 'LineWidth', 2, 'MarkerSize', 10, 'DisplayName', 'goal position');
plot(x(1,:), x(2,:), 'g-', 'LineWidth', 2, 'DisplayName', 'optimal trajectory');
plot(x_s(:,1), x_s(:,2), 'r--', 'LineWidth', 2, 'DisplayName', 'simulated trajectory');

hold off;
xlabel('x');
ylabel('y');
title(['vehicle trajectory through ' current_model ' current']);
legend();
axis('square');

set(findall(gcf, '-property', 'FontSize'), 'FontSize', 16);
set(findall(gcf, '-property', 'FontWeight'), 'FontWeight', 'bold');

%% compute cost function to minimize control input
function [J] = compute_cost(z, N)
[~, u, tauf] = extract_state(z, N);

w = lpm.compute_quadrature_weights(N);
L = 0.5 * sum(u.^2, 1);
J = (tauf / 2) * sum(w .* L);
end

%% initialize decision variables
function [z0] = initialize_guess(N, x0, xf)
% final time guess
tauf = 100;

% heading toward goal (atan2 of east/north for compass convention)
theta0 = atan2(xf(1) - x0(1), xf(2) - x0(2));

% straight line interpolation for position, constant heading
tau = lpm.compute_collocation_points(N);
s   = 0.5 * (tau + 1);
xy  = x0(1:2) * (1 - s) + xf * s;
x   = [xy; theta0 * ones(1, N+1)];

% constant speed, zero steering
u = [ones(1, N+1); zeros(1, N+1)];

z0 = [x(:); u(:); tauf];
end

%% extract state variables from decision vector
function [x, u, tauf] = extract_state(z, N)
ix = 3 * (N + 1);
iu = ix + 2 * (N + 1);

x    = reshape(z(1:ix), 3, []);
u    = reshape(z(ix+1:iu), 2, []);
tauf = z(end);
end

%% compute constraints
function [c, ceq] = nonlcon(z, x0, xf, current_field, wheelbase, N)
% extract state variables
[x, u, tauf] = extract_state(z, N);

% compute collocation points
tau = 0.5 * tauf * (lpm.compute_collocation_points(N) + 1);

% compute differentiation matrix
D = lpm.compute_differentiation_matrix(N);

% preallocate output
c   = [];
ceq = zeros(3 * (N + 1) + 5, 1);

% compute equality constraints due to system dynamics
ic        = 3 * (N + 1);
ceq(1:ic) = reshape(x * D' - 0.5 * tauf * system_dynamics(tau, x, u, current_field, wheelbase), [], 1);

% compute equality constraints due to boundary conditions
ceq(ic+1:ic+3) = x(:,1) - x0;
ceq(ic+4:ic+5) = x(1:2,end) - xf;
end

%% system dynamics of kinematic vehicle with current
function [xdot] = system_dynamics(~, x, u, current_field, wheelbase)
theta = x(3,:);
v     = u(1,:);
delta = u(2,:);
c     = current_field(x(1,:), x(2,:));

xdot = [v .* sin(theta) + c(1,:);
        v .* cos(theta) + c(2,:);
        (v / wheelbase) .* tan(delta)];
end