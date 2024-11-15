% Implements an unscented Kalman Filter (goole "Sigma-point Kalman Filters")
% This code is derived from: https://stackoverflow.com/questions/55813719/multi-sensors-fusion-using-kalman-filter

function [] = ukf()
    % time step
    dt = 0.01;

    t=(0:dt:2)';
    n = numel(t);

    %ground truth
    signal = sin(t)+t;

    % state matrix: (pos, vel, acc)
    X = zeros(3,1);

    % covariance matrix
    P = zeros(3,3);

    % kalman filter output through the whole time
    X_arr = zeros(n, 3);

    % system noise
    Q = [0.04 0 0;
         0 0 0;
         0 0 0];

    % transition matrix, a way to represent the process, and control. Only
    % works for linear processes, will be replaced soon.
    F = [1 dt 0;
         0 1 0;
         0 0 1];

    % observation matrix
    H = [1 0 0];

    % variance of signal 1
    s1_var = 0.08*ones(size(t));
    s1 = generate_signal(signal, s1_var);

    % variance of signal 2
    s2_var = 0.01*(cos(8*t)+10*t);
    s2 = generate_signal(signal, s2_var);

    % variance of signal 3
    s3_var = 0.02*(sin(2*t)+2);
    s3 = generate_signal(signal, s3_var);

    % variance of signal 4
    s4_var = 0.06*ones(size(t));
    s4 = generate_signal(signal, s4_var);

    % fusion
    for i = 1:n
        if (i == 1)
            [X, P] = init_kalman(X, s1(i, 1)); % initialize the state using the 1st sensor
        else
            % Prediction stage
            [xx, ww] = calculate_sigma_points_vdm(X, P);
            % Propagate sigma-points
            for i = 1:size(xx)(1,2)
                xx(:,i) = prediction(xx(:,i), P, Q, F);
            end
            [X, P] = unscented_transform(xx, ww, Q);

            [X, P] = unscented_update(xx, ww, X, P, H, s1(i, 1), [s1(i, 2) 0 0]);
            % [X, P] = update(X, P, s1(i, 1), s1(i, 2), H);
            % [X, P] = update(X, P, s2(i, 1), s2(i, 2), H);
            % [X, P] = update(X, P, s3(i, 1), s3(i, 2), H);
            % [X, P] = update(X, P, s4(i, 1), s4(i, 2), H);
        end

        X_arr(i, :) = X;
    end

    plot(t, signal, 'LineWidth', 4);
    hold on;
    plot(t, s1(:, 1), '--', 'LineWidth', 1);
    plot(t, s2(:, 1), '--', 'LineWidth', 1);
    plot(t, s3(:, 1), '--', 'LineWidth', 1);
    plot(t, s4(:, 1), '--', 'LineWidth', 1);
    plot(t, X_arr(:, 1), 'LineWidth', 4);
    hold off;
    grid on;
    legend('Ground Truth', 'Sensor Input 1', 'Sensor Input 2', 'Sensor Input 3', 'Sensor Input 4', 'Fused Output');
end

function [X, P] = unscented_update(yy, ww, X, P, H, Z, R)
    % - H: observation matrix. Just a way to represent conversions from
    %   measurements to actual data, will become obsolete
    % Read: https://docs.duckietown.com/daffy/course-intro-to-drones/ukf/theory/ukf-specifics.html
    % - yy: predicted sigma points after unscented transform
    % - R: noise covariance matrix
    % - Z: measurement matrix (signal from the sensor)
    % - X: predicted state after unscented transform
    npoints = size(yy)(1, 2);
    dim = size(yy)(1, 1);

    % Compute the measurement sigma points
    zz = zeros(size(yy));
    for i = 1:npoints
        % TODO: Matrix H will become obsolete, it is not applicable to a general
        % case, because the observation coversion cannot be expressed as a
        % linear transformation at all times
        zz(:, i) = yy(:, i);
    end
    % Compute mean of the measurement sigma points using unscented transform
    mu = zeros(dim, 1);
    for i = 1:npoints
        mu = mu + ww(i) * zz(:, i);
    end
    % Compute covariance of the measurement sigma points using unscented transform
    Pz = zeros(dim, dim);
    for i = 1:npoints
        Pz = ww(i) * (zz(:, i) - mu) * (zz(:, i) - mu)';
    end
    Pz = Pz + R;
    Pz, inv(Pz)

    % Compute residual in measurement space
    yt = Z - mu;
    % Compute cross-covariance between state and measurements
    Pxz = zeros(dim, dim);
    for i = 1:npoints;
        Pxz = Pxz + ww(i) * (yy(:, i) - X) * (zz(:, i) - mu)';
    end
    % Compute Kalman gain
    K = Pxz * inv(Pz);

    % Compute the mean of posterior state estimate
    X = X + K * yt;
    % Compute covariance matrix of posterior state estimate
    P = P - K * Pz * K';
end

function [X, P] = unscented_transform(yy, ww, Q)
    % Ins:
    % - xx: propagated sigma-points (post-prediction), representation is the
    % same as the one used for xx in `calculate_sigma_points_vdm`
    % - ww: weights, see `calculate_sigma_points_vdm`
    % - Q: system noise
    %
    % Outs:
    % - X: prior mean
    % - P: covariance
    %
    % See https://docs.duckietown.com/daffy/course-intro-to-drones/ukf/theory/ukf-specifics.html

    % Calculate prior mean
    npoints = size(yy)(1, 2); % number of sigma-points
    dim = size(yy)(1, 1);
    X = zeros(dim, 1);
    for i = 1:npoints
        X = X + ww(i) * yy(:,i);
    end

    % Calculate covariance matrix
    P = zeros(dim, dim);
    for i = 1:npoints
        P = P + ww(i) * (yy(:, i) - X) * (yy(:,i) - X)';
    end
    P = P + Q;
end

function [s] = generate_signal(signal, var)
    noise = randn(size(signal)).*sqrt(var);

    s(:, 1) = signal + noise;

    % Just a way to store square of variance, it has no special meaning
    s(:, 2) = var;
end

function [X, P] = init_kalman(X, y)
    X(1,1) = y;
    X(2,1) = 0;

    P = [100 0 0;
         0   300 0;
        0 0 0];
end

function [xx, ww] = calculate_sigma_points_vdm(X, P)
    % VDM stands for "Van der Merwe"
    % Calculate sigma points using Van der Merwe's et. al. algorithm (2004)
    % https://www.researchgate.net/publication/228846510_Sigma-Point_Kalman_Filters_for_Nonlinear_Estimation_and_Sensor-Fusion-Applications_to_Integrated_Navigation
    % xx - a matrix whose columns represent sigma functions, and
    % ww - weights, same structure as xx
    % zeta - scaling factor
    % The number of sigma points N = 2 * dim + 1, where dim - number of dimensions of the state vector

    dim = size(P)(1, 1);
    % Inverse square root of a matrix
    psqr = sqrtm(P);
    % TODO: I have no idea what zeta should be equal to
    zeta = 0.1;

    % Initialize sigma-points
    xx = zeros(dim, 2 * dim + 1);
    xx(:,1) = X;
    for i = 1:dim
        xx(:, i + 1) = X + zeta * psqr(:, i);
    end
    for i = 1:dim
        xx(:, i + dim + 1) = X - zeta * psqr(:, i);
    end
    % Initialize weights. Results in normalized weight vector
    ww = zeros(dim, 2 * dim + 1);
    ww = [1, ones(1, dim * 2)];
    ww = ww / sum(ww);
end

function [X, P] = prediction(X, P, Q, F)
    X = F*X;
    P = F*P*F' + Q;
end

function [X, P] = update(X, P, y, R, H)
    % R: measurement noise covariance matrix (in 1-dimensional case, may just be represented as a scalar)
    % Difference b/w measured, and previous, in measurement space
    Inn = y - H*X;

    S = H*P*H' + R;
    K = P*H'/S;

    X = X + K*Inn;
    P = P - K*H*P;
end
