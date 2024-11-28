addpath("ukf")
pkg load matgeom;
source("util/util.m")
source("imu_estimate.m")

[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

% Imu calculations
% Z vector in global frame. The UAV is stable, acc. measurements are aligned with the gravitational vector (G-vector)
gvec = [ttyyaccgyro(1, 3), ttyyaccgyro(1, 4), ttyyaccgyro(1, 5)];
rot = imu_rotation_from_acc(gvec);
grid on
hold on
%plot3_rot(rot)
euler = rot2euler(rot);
yaw = 0;
hold on;
load yy.msession;

function [alt] = predict_baro(alt0, uarg)
	% Estimate speed

	dt = uarg.xx(2, 1) - uarg.xx(1, 1);
	if uarg.i >= 3
		dx = uarg.xx(uarg.i - 1, 3) - uarg.xprev;
		v = dx / dt;
	else
		v = 0;
	end
	% Accelerometer value
	persistent accpos = 0; % Position hint for acc measurements
	t = uarg.xx(uarg.i, 1);
	apos = ts_find_measurement_before(uarg.acct, t);
	at = uarg.acct(apos);
	az = uarg.accz(apos);
	alt = predict_altitude_taylor(uarg.xprev, dt, v, az);
end

% SR-UKF using only barometric data
if true
	% Implement SR-UKF
	n = 1; % Dimension of state
	m = 1; % Dimension of measurement
	q = 0.1; % Std of process
	r = 0.1; % Std of measurement
	Qs = q * eye(n); % std matrix of process
	Rs = r * eye(m); % std matrix of measurement
	h = @(x, ~)[x]; % TODO: measurement function
	f = @(x, a)predict_baro(x, a);
	x = ttyygnss(1, 3); % Initial state
	S = eye(n); % initial square root of state covariance
	N = size(ttyygnss)(1); % total dynamic steps
	N
	%N = 10000; % total dynamic steps
	xV = zeros(n, N); % estmate
	barowinsize = 10; % Use "barowinsize" last measurements to estimate standard deviation
	for k = 2:N
		% Estimate baro variance by 10 last available measurements
		Rs = std(ttyybaroalt(max(1, k - barowinsize):k));
		% Estimate process varaince by 10 last assessments (WARNING: won't scale for d>1 dimensions, diagonal of Qs must be set)
		Qs = std(xV(1, max(1, k - barowinsize):k));
		z = ttyybaroalt(k, 3); % TODO: get the next measurement
		context.xx = ttyybaroalt;
		context.i = k;
		context.xprev = xV(:, k - 1);
		% Has been calculated previously, saved in a file
		context.accz = yy;
		context.acct = tt;
		[x, S] = srukf(f, x, S, h, z, Qs, Rs, context); % ukf
		xV(:, k) = x; % save estimate
		k
	end
	% TODO: plot
end

if false
	load yy
	% Plot data
	% Baro
	grid on
	ax1 = subplot(2, 1, 1)
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3);
	plot(ttbaro, yybaro)
	% Accel
	ax2 = subplot(2, 1, 2)
	plot(tt, yy)
	linkaxes([ax1 ax2], 'x')

	% FFT
	% Sampling freq.
	f = 1 / (tt(2) - tt(1));
	% period
	t = 1 / f;
	% length
	l = numel(yy);

	yfft = fft(yy);
	% plot
	figure
	plot(f / l * (0:l - 1), abs(yfft));
	title("Complex Magnitude of fft Spectrum")
	xlabel("f (Hz)")
	ylabel("|fft(X)|")
end

if false % Plot vertial acceleration graph along w/ baro alt
	tt = zeros(size(ttyyaccgyro)(1), 1);
	yy = zeros(size(ttyyaccgyro)(1), 1);
	t = ttyyaccgyro(1, 1);
	tt(1) = t;
	for i = 2:size(ttyyaccgyro)(1)
		% Prepare data
		acc = ttyyaccgyro(i, 3:5);
		gyro = ttyyaccgyro(i, 7:9);
		prevt = t;
		t = ttyyaccgyro(i, 1);
		dt = t - prevt;
		% Update rotation matrix in global frame
		euler = euler_estimate_complemetary(euler, acc, gyro, dt);
		i / size(ttyyaccgyro)(1)
		% Get gravitational component
		az = vert_accel(euler, acc');
		% Fix data
		tt(i) = t;
		yy(i) = az;
	end

	% Plot data
	% Baro
	grid on
	ax1 = subplot(2, 1, 1)
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3);
	plot(ttbaro, yybaro)
	% Accel
	ax2 = subplot(2, 1, 2)
	plot(tt, yy)
	linkaxes([ax1 ax2], 'x')
end

if false % Animated plot of position estimation
	startpos = 1
	printeach = 10
	updateeach = printeach * 5
	t = ttyyaccgyro(startpos, 1)
	for i = startpos:size(ttyyaccgyro)(1)
		acc = ttyyaccgyro(i, 3:5);
		gyro = ttyyaccgyro(i, 7:9);
		prevt = t;
		t = ttyyaccgyro(i, 1);
		dt = t - prevt;
		euler = euler_estimate_complemetary(euler, acc, gyro, dt);
		if mod(i, updateeach) == 0
			cla
		end
		if mod(i, printeach) == 0
			hold on;
			xlim([-0.3 0.3])
			ylim([-0.3 0.3])
			zlim([-0.3 0.3])
			plot3_euler(euler);
			pause(0.02)
			t
		end
	end
end

if false % Plot raw data
	hold on;
	grid on;
	% Baro
	ax1 = subplot(2, 1, 1);
	ttbaro = ttyybaroalt(:, 1)
	yybaro = ttyybaroalt(:, 3)
	plot(ttbaro, yybaro);
	% GNSS
	ax2 = subplot(2, 1, 2);
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);
	plot(ttgnss, yygnss)
	linkaxes([ax1, ax2], "x")
	% Acc / Gyro
	figure
	subplot(2, 1, 1)
	hold on
	grid on
	ttacc = ttyyaccgyro(:, 1);
	for i = [3 4 5]
		plot(ttacc, ttyyaccgyro(:, i));
	end
	subplot(2, 1, 2)
	for i = [7 8 9]
		hold on
		grid on
		plot(ttacc, ttyyaccgyro(:, i));
	end
end
