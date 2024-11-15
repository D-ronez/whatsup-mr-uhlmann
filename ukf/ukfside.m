function ukf_test()
	% Initial states
	y = [0; 0];
	u = [1; 2];
	xhat = [0; 0];
	P = [5  0; 0 2];
	Q = [1, 0; 0, 2];
	R = [1.5, 0; 0, 2];
	a = 1;
	k = 2;
	b = 3;
	M = 2;
	for i = 1:2
		[xhat, y, P] = ukf(xhat, y, u, P, Q, R, a, k, b, M);
	end
end

function [xhat, y, P] = ukf(xhat, y, u, P, Q, R, a, k, b, M)
	column = 2 * M + 1;
	row = M;

	% Step 0 - Create the weights
	[WM, Wc] = ukf_create_weights(a, b, k, row);

	% UPDATE: Step 1 - Compute sigma points
	[xhati] = ukf_compute_sigma_points(xhat, P, a, k, row);

	% UPDATE: Step 2 - Use the nonlinear measurement function to compute the predicted measurements for each of the sigma points.
	yhati = xhati; % Here we assume that the observation function y = h(x, u) = x

	% UPDATE: Step 3 - Combine the predicted measurements to obtain the predicted measurement
	yhat = ukf_multiply_weights(yhati, WM, M);

	% UPDATE: Step 4 - Estimate the covariance of the predicted measurement
	Py = ukf_estimate_covariance(yhati, yhat, Wc, R, row);

	% UPDATE: Step 5 - Estimate the cross-covariance between xhat and yhat. Here i begins at 1 because xhati(0) - xhat(0) = 0
	Pxy = ukf_estimate_cross_covariance(xhati, xhat, yhati, yhat, a, k, M);

	% UPDATE: Step 6 - Find kalman K matrix
	K = ukf_create_kalman_K(Py, Pxy, M);

	% UPDATE: Step 7 - Obtain the estimated state and state estimation error covariance at time step
	[xhat, P] = ukf_state_update(K, Py, P, xhat, y, yhat, M);

	% PREDICT: Step 0 - Predict the state and state estimation error covariance at the next time step
	xhati = ukf_compute_sigma_points(xhat, P, a, k, M);

	% PREDICT: Step 1 - Use the nonlinear state transition function to compute the predicted states for each of the sigma points.
	xhati = ukf_transition(xhati, u, M);

	% PREDICT: Step 2 - Combine the predicted states to obtain the predicted states at time
	xhat = ukf_multiply_weights(xhati, WM, M);

	% PREDICT: Step 3 - Compute the covariance of the predicted state
	P = ukf_estimate_covariance(xhati, xhat, Wc, Q, M);

end

function [WM, Wc] = ukf_create_weights(a, b, k, M)
	column = 2 * M + 1;
	WM = zeros(1, column);
	Wc = zeros(1, column);
	for i = 1:column
		if(i == 1)
			Wc(i) = (2 - a * a + b) - M / (a * a * (M + k));
			WM(i) = 1 - M / (a * a * (M + k));
		else
			Wc(i) = 1 / (2 * a * a * (M + k));
			WM(i) = 1 / (2 * a * a * (M + k));
		end
	end
end

function [xi] = ukf_compute_sigma_points(x, P, a, k, M)
	column = 2 * M + 1;
	compensate_column = 2 * M - 1;
	row = M;
	c = a * a * (M + k);
	xi = zeros(row, column);

	% According to the paper "A New Extension of the Kalman Filter to Nonlinear Systems"
	% by Simon J. Julier and Jeffrey K. Uhlmann, they used L = chol(c*P) as "square root",
	% instead of computing the square root of c*P. According to them, cholesky decomposition
	% was a numerically efficient and a stable method.
	L = chol(c*P, 'lower');

	for j = 1:column
		if(j == 1)
			xi(:, j) = x;
		elseif(and(j >= 2, j <= row + 1))
			xi(:, j) = x + L(:, j - 1);
		else
			xi(:, j) = x - L(:, j - compensate_column);
		end
	end

end

function x = ukf_multiply_weights(xi, W, M)
	column = 2 * M + 1;
	row = M;
	x = zeros(row, 1);
	for i = 1:column
		x = x + W(i)*xi(:, i);
	end
end

function P = ukf_estimate_covariance(xi, x, W, O, M)
	column = 2 * M + 1;
	row = M;
	P = zeros(row, row);
	for i = 1:column
		P = P + W(i)*(xi(:, i) - x)*(xi(:, i) - x)';
	end
	P = P + O;
end

function P = ukf_estimate_cross_covariance(xi, x, yi, y, a, k, M)
	column = 2 * M + 1;
	row = M;
	c = 1 / (2 * a * a * (M + k));
	P = zeros(row, row);
	for i = 2:column % Begins at 2 because xi(:, 1) - x = 0
		P = P + (xi(:, i) - x)*(yi(:, i) - y)';
	end
	P = c*P;
end

function K = ukf_create_kalman_K(Py, Pxy, M)
	row = M;
	K = zeros(row, row);
	for i = 1:row
		% Solve Ax = b with Cholesky
		L = chol(Py, 'lower');
		y = linsolve(L, Pxy(:, i));
		K(:, i) = linsolve(L', y);
	end
	% This will work to K = linsolve(Py, Pxy);
end

function [xhat, P] = ukf_state_update(K, Py, P, xhat, y, yhat, M)
	row = M;
	xhat = xhat + K*(y - yhat);
	P = P - K*Py*K';
end

function xhati = ukf_transition(x, u, M)
	column = 2 * M + 1;
	row = M;
	xhati = zeros(row, column);
	for i = 1:column
		xhati(:, i) = transistion_function(x(:, i), u);
	end
end

function dx = transistion_function(x, u)
	dx = zeros(2, 1);
	dx(1) = -2*x(1)*x(2) + 4*x(2) + 4*u(1);
	dx(2) = -x(1) - 3*x(2) + 7*u(2);
end
