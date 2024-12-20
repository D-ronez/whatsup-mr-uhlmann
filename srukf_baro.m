function [] = srukf_baro()
	addpath("ukf");
	% IMU assessment, projection onto Z axis is already computed, load the values
	load tt.msession;
	load az.msession;
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");
	r_baro = 1.2; % std of baro measurement
	r_gnss = 0.1 % std of GNSS measurement
	q_process = 0.1; % std of the process
	% Prediction function
	f = @(x, uarg)[predict_altitude_taylor(x, uarg.dt, uarg.av, uarg.az)];
	% Measurement function
	h = @(x, ~)[x];
	% Sensor measurements, baro, accel's Z post IMU estimation
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3);
	ttaz = tt;
	yyaz = az;
	% Covariance matrix
	s = 1;

	% Compensate for freefall acceleration
	g = 9.80665;

	yyfuse = zeros(size(yybaro));
	n = size(ttbaro)(1);
	barooffset = ttyygnss(1, 5)
	x = yybaro(1) + barooffset;
	uarg.dt = ttbaro(2) - ttbaro(1);
	posaa = 1;
	posgnss = 1;
	yyfuse(1) = x;
	for k = 2:n
		% Prepare arg-s for the prediction function.
		t = ttbaro(k);
		posaa = get_tt_pos(ttbaro, t, posaa);
		uarg.dt = ttbaro(k) - ttbaro(k - 1);
		uarg.az = yyaz(posaa) + g;
		uarg.av = (yybaro(k) - yybaro(k - 1)) / (uarg.dt);

		posgnss = get_tt_pos(ttyygnss(:, 1), t, posgnss);
		usegnss = true;
		if isnan(posgnss)
			usegnss = false;
			posgnss = 1;
		end
		altgnss = ttyygnss(posgnss, 5);

		altbaro = yybaro(k) + barooffset;

		% Run SR-UKF
		z = [altbaro, altgnss];
		r = [r_baro, r_gnss];
		if ~usegnss
			z = z(1)
			r = r(1)
		end
		% [x, s] = srukf(f, x, s, h, z, q_process, r_baro, uarg);
		[x, s] = ukf1d(f, h, s, r, q_process, x, z, uarg);

		% Save the result
		yyfuse(k) = x;
		k
	end

	save res;

	plot(ttbaro, yyfuse);
	hold on
	plot(ttbaro, yybaro + barooffset);
	plot(ttyygnss(:, 1), ttyygnss(:, 5));
	legend('fuse', 'baro', 'gnss')
end

function [pos] = get_tt_pos(tt, t, hint)
	% Gets appropriate position for timeseries w/ different timestep
	tt = tt(hint:end, 1);
	pos = nan;
	for i = 1:size(tt)(1)
		if t < tt(i)
			pos = i + hint;
			break
		end
	end
end

function [alt1] = predict_altitude_taylor(alt0, dt, av0, az0)
	% Part of Kalman filter. Predicts altitude using Taylor series
	% "Midpoint Method".
	%
	% alt0: Z-altitude at previous time
	% dt:   time delta between IMU measurements
	% az0:  [m/s^2] Accelerometer's Z reading
	% alt1: predicted altitude

	alt1 = alt0 + dt * av0 + dt^2 * az0 / 2;
end
