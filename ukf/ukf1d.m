function [x, p] = ukf1d(f, h, p, r, xprev, z, uarg)
	% 1-dimensional version of Unscented Kalman filter
	% f - prediction function @(xprev, uarg)
	% h - measurement function @(xprev, uarg)
	% p - variance of process
	% r - variance of measurement
	% z - measurement
	%

	m = numel(z); %numer of measurements
	alpha = 1e-3; %default, tunable
	ki = 0; %default, tunable
	beta = 2; %default, tunable
	l = numel(xprev);
	lambda = alpha^2 * (l + ki) - l; %scaling factor
	c = l + lambda; %scaling factor

	% weights
	wm = [lambda / c 0.5 / c + zeros(1, 2 * l)]; %weights for means
	wc = wm;
	wc(1) = wc(1) + (1 - alpha^2 + beta); %weights for covariance

	% Apriori
	xx = sigmas(xprev, p, c);
	% Predict
	for i = 1:size(xx, 2)
		xx(:, i) = f(xx(:, i), uarg);
	end
	xk = sum(xx .* wm); % Predicted x (priori)
	p = sum(wc .* [xx - xk].^2) + r;

	% Posteriori
	yy = xx;
	for i = 1:size(xx, 2)
		yy(:, i) = h(yy(:, i), uarg);
	end
	yk = sum(wm .* yy);

	p_ykyk = sum(wc .* [yy - yk].^2) + r;
	p_xkyk = sum(wc .* [xx - xk] .* [yy - yk]);
	k = p_xkyk / p_ykyk; % Kalman gain
	x = xk + k * (z - yk);
	p = p - k * p_ykyk * k;
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
