function [] = srukf_baro()
	addpath("ukf");
	% IMU assessment, projection onto Z axis is already computed, load the values
	load tt.msession;
	load az.msession;
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");
	% STD values
	r_gnss = 0.005;
	q_process = 12.0;
	% Correction scale for baro offset, the lower it is, the slower baro values are corrected
	bcorr_slowdown = 0.008;
	% Prediction function
	f = @(x, uarg)[predict_altitude_taylor(x, uarg.dt, uarg.av, uarg.az)];
	% Measurement function
	h = @(x, ~)[x];
	% Baro altitude timeseries
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3);
	% Gnss altitude timeseries
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);
	ttaz = tt;
	yyaz = az;
	% Covariance matrix
	s = r_gnss;

	% Compensate for freefall acceleration
	g = 9.80665;
	azcalib = 2.35;
	%g = 0;

	n = numel(ttgnss);
	x = yygnss(1);
	tprev = ttgnss(1);
	t = tprev
	uarg.dt = 0;
	posaccel = 1;
	posbaro = 1;
	posgnss = 1;
	yyfuse = [x];
	ttfuse = [t];
	ttmse = [];
	yymse = [];
	ttmae = [];
	yymae = [];
	ttbarocorr = [];
	yybarocorr = [];

	% Compensate baro drift using GNSS
	barooffset = ttyygnss(1, 5);
	ib = 1;
	ig = 1;
	i = 1;
	yb = 0;
	yg = 0;
	while ~isnan(ib) && ~isnan(ig)
		% Previous sensor measurement
		yb = yybaro(ib) + barooffset;
		yg = yygnss(ig);

		% Model acquisition of the next value from the pre-loaded time series
		[ig, ib, used] = pick_next_point_t(ig, ib, ttgnss, ttbaro);
		usegnss = (used == 1);
		t = 0;
		if isnan(ib) || isnan(ig)
			break
		end

		% Next sensor measurement
		if usegnss
			gnssready = true;
			yg = yygnss(ig);
			t = ttgnss(ig);

			% Update baro offset, apply correction with each GNSS sensor update
			barooffsetdelta = (yg - yb) * bcorr_slowdown;
			barooffset = barooffset + barooffsetdelta;
		else
			baroready = true;
			yb = yybaro(ib);
			t = ttbaro(ib);

			% Save corrected barometer values
			yb = yybaro(ib) + barooffset;
			ttbarocorr = [ttbarocorr, t];
			yybarocorr = [yybarocorr, yb];
			yybaro(ib - 1) = yb;
		end

		i = i + 1
	end

	for k = 2:n
		t = ttgnss(k)
		ygnss = yygnss(k);

		% Prepare arg-s for the prediction function.
		posaccel = get_tt_pos(ttaz, tprev, posaccel) - 1
		assert(posaccel > 0);
		winstart = max([1, posaccel - 0]);
		taz = ttaz(posaccel)
		%uarg.az = yyaz(posaccel) + g + azcalib;
		uarg.az = 0;
		uarg.dt = t - tprev;

		% Calculate speed using linear approximation
		posbaro = get_tt_pos(ttbaro, tprev, posbaro);
		winstart = max([1, posbaro - 50]);
		sparse_factor = 5;
		tt = ttbaro(winstart:posbaro);
		yy = yybaro(winstart:posbaro);
		uarg.av = polyfit(tt(1:sparse_factor:numel(tt)), yy(1:sparse_factor:numel(yy)), 1)(1)
		% Baro start
		if false
			posbaro = get_tt_pos(ttbaro, tprev, posbaro)
			tbaro = ttbaro(posbaro)
			if posbaro > 1
				dybaro = yybaro(posbaro) - yybaro(posbaro - 1)
				dtbaro = ttbaro(posbaro) - ttbaro(posbaro - 1)
				% Calculate speed using baro
				uarg.av = dybaro / dtbaro
			else
				uarg.av = 0
			end
		end

		% Run SR-UKF
		z = ygnss;
		r = r_gnss;
		[x, s] = ukf1d(f, h, s, r, q_process, x, z, uarg);

		% Save the result
		yyfuse(k) = x;
		ttfuse(k) = t;

		% Save MSE
		msewinsize = 100;
		winend = numel(yyfuse);
		winbegin = max(1, winend - msewinsize);
		a = yyfuse(winbegin:winend);
		b = yygnss(winbegin:winend)';
		ttmse = [ttmse, t];
		yymse = [yymse, sum(a - b) .^ 2 / msewinsize];
		% Save MAE
		ttmae = [ttmae, t];
		yymae = [yymae, sum(abs(a - b))];

		tprev = t;
		k
	end

	save ukfgnss;

	figure;
	hold on
	plot(ttbaro, yybaro);
	plot(ttgnss, yygnss);
	plot(ttfuse, yyfuse);
	legend('baro', 'gnss', 'fuse');

	% Plot MSE
	figure;
	plot(ttmse, yymse);
end

function [pos] = get_tt_pos(tt, t, hint)
	% Gets appropriate position for timeseries w/ different timestep
	tt = tt(hint:end, 1);
	pos = nan;
	for i = 1:size(tt)(1)
		if t < tt(i)
			pos = i + hint - 1;
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

function [ia, ib, i] = pick_next_point_t(ia, ib, tta, ttb)
	% From 2 sets of time points, pick up the next one (closer to the current)
	%

	i = 1;
	if ia >= numel(tta)
		ia = nan;
	elseif ib >= numel(ttb)
		ib = nan;
	else
		ta = tta(ia + 1);
		tb = ttb(ib + 1);
		if ta > tb
			ib = ib + 1;
			i = 2;
		else
			ia = ia + 1;
			i = 1;
		end
	end
end
