% clear
% clc
% close all

dt = 0.1;
T = 120;
N = T/dt;

%% State vector
% [x y z vx vy vz bax bay baz]
n = 9;

%% True values
pos = [0;0;0];
vel = [2;0.5;0.2];

true_acc_bias = [0.02;-0.015;0.01];

%% EKF initial estimate
x_est = zeros(n,1);

P = eye(n)*1;

%% Noise models
imu_noise = 0.02;
gps_noise = 1.5;

Q = diag([0 0 0 0.01 0.01 0.01 1e-6 1e-6 1e-6]);
R = eye(3)*gps_noise^2;

%% Visualization
figure
grid on
hold on
axis equal
view(3)
xlabel('X')
ylabel('Y')
zlabel('Z')

true_traj = [];
est_traj = [];

for k = 1:N

    %% TRUE MOTION
    acc_true = [0;0;0];

    vel = vel + acc_true*dt;
    pos = pos + vel*dt;

    %% IMU measurement
    imu_acc = acc_true + true_acc_bias + imu_noise*randn(3,1);

    %% GPS availability
    gps_available = (k < 400) || (k > 700);

    if gps_available
        gps_meas = pos + gps_noise*randn(3,1);
    end

    %% EKF PREDICTION

    x = x_est;

    pos_est = x(1:3);
    vel_est = x(4:6);
    bias_est = x(7:9);

    acc_est = imu_acc - bias_est;

    pos_pred = pos_est + vel_est*dt + 0.5*acc_est*dt^2;
    vel_pred = vel_est + acc_est*dt;

    x_pred = [pos_pred; vel_pred; bias_est];

    %% Jacobian
    F = eye(n);

    F(1:3,4:6) = eye(3)*dt;
    F(1:3,7:9) = -0.5*eye(3)*dt^2;
    F(4:6,7:9) = -eye(3)*dt;

    P_pred = F*P*F' + Q;

    %% MEASUREMENT UPDATE (GPS)

    if gps_available

        H = zeros(3,n);
        H(:,1:3) = eye(3);

        K = P_pred*H'/(H*P_pred*H' + R);

        x_est = x_pred + K*(gps_meas - H*x_pred);

        P = (eye(n)-K*H)*P_pred;

    else

        x_est = x_pred;
        P = P_pred;

    end

    %% Save trajectories
    true_traj = [true_traj pos];
    est_traj = [est_traj x_est(1:3)];

    %% Plot
    clf
    hold on
    grid on
    axis equal
    view(3)

    plot3(true_traj(1,:),true_traj(2,:),true_traj(3,:),'b','LineWidth',2)
    plot3(est_traj(1,:),est_traj(2,:),est_traj(3,:),'r--','LineWidth',2)

    %% Covariance ellipsoid
    P_pos = P(1:3,1:3);
    plot_ellipsoid(x_est(1:3),P_pos,2)

    if gps_available
        title('GPS Available')
    else
        title('GPS LOST (INS drifting)')
    end

    drawnow

end


%% Ellipsoid plotting function
function plot_ellipsoid(mu,P,k)

[V,D] = eig(P);

[xs,ys,zs] = sphere(20);
pts = [xs(:) ys(:) zs(:)]';

ell = V*sqrt(D)*pts*k;

ell(1,:) = ell(1,:) + mu(1);
ell(2,:) = ell(2,:) + mu(2);
ell(3,:) = ell(3,:) + mu(3);

X = reshape(ell(1,:),size(xs));
Y = reshape(ell(2,:),size(ys));
Z = reshape(ell(3,:),size(zs));

surf(X,Y,Z,'FaceAlpha',0.3,'EdgeColor','none')

end