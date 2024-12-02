function [x, S] = srukfcustom(fstate, x, S, hmeas, z, Qs, Rs, uarg)
	% multi-dimensional version of UnscentedKalman filter
	%
	% SR-UKF   Square Root Unscented Kalman Filter for nonlinear dynamic systems
	% [x, S] = ukf(f, x, S, h, z, Qs, Rs) returns state estimate, x and state covariance, P
	% for nonlinear dynamic system (for simplicity, noises are assumed as additive):
	%           x_k + 1 = f(x_k) + w_k
	%           z_k   = h(x_k) + v_k
	% where w ~ N(0, Q) meaning w is gaussian noise with covariance Q
	%       v ~ N(0, R) meaning v is gaussian noise with covariance R
	% Inputs:   f: function handle for f(x)
	%           x: "a priori" state estimate
	%           S: "a priori" estimated the square root of state covariance
	%           h: fanction handle for h(x)
	%           z: current measurement
	%           Qs: process noise standard deviation
	%           Rs: measurement noise standard deviation
	%           hmeas: measurement function
	%           fstate: prediction function
	%           uarg: user argument, context that will be passed to measurement, and prediction function
	% Output:   x: "a posteriori" state estimate
	%           S: "a posteriori" square root of state covariance
	%
	% Example:
	% Reference: R. van der Merwe and E. Wan.
	% The Square-Root Unscented Kalman Filter for State and Parameter-Estimation, 2001
	%
	% By Zhe Hu at City University of Hong Kong, 05 / 01 / 2017
	%

	L = numel(x); %numer of states

	% Scaling parameters
	m = numel(z); %numer of measurements
	alpha = 1e-3; %default, tunable
	ki = 0; %default, tunable
	beta = 2; %default, tunable
	lambda = alpha^2 * (L + ki) - L; %scaling factor
	c = L + lambda; %scaling factor

	% Weights
	Wm = [lambda / c 0.5 / c + zeros(1, 2 * L)]; %weights for means
	Wc = Wm;
	Wc(1) = Wc(1) + (1 - alpha^2 + beta); %weights for covariance

	% Sigmas
	c = sqrt(c);
	X = sigmas(x, S, c); %sigma points around x
	for i = 1:size(X, 2):
		X(:, i) = fstate(X(:, i), uarg)
	end

	% Sigma point calculation, and time update
end

function X = sigmas(x, S, c)
	%Sigma points around reference point
	%Inputs:
	%       x: reference point
	%       S: square root of covariance
	%       c: coefficient
	%Output:
	%       X: Sigma points

	A = c * S';
	Y = x(:, ones(1, numel(x)));
	X = [x, Y + A, Y-A];
end
