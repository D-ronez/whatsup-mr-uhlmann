% Acquires fused altitude value from baro, and GNSS, and performs scaled
% on-the-fly correction of baro using GNSS to prevent it from drifting away over
% time.
%
% Variance-based filter for generating fused altitude from GNSS, and baro
% measurements. It computes trust scores for GNSS, and baro, and produces fused
% value based on weighted sum of both sensors whose weights are inferred
% according to their respective trust values.
%

function [x, istate] = coolfilter(z, i, istate)
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

	barooffset = 0; % The baro values will be constantly corrected against GNSS
	bcorr_slowdown = 0.008; % Correction scale for baro offset, the lower it is, the slower baro values are corrected
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3) + barooffset;;
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);

	% Cooldown for accumulated error
	intcooldown = 0.2; % [m/s]
	bsigma = 0.9;
	gsigma = 0.1;
	windowsize = 25;

	ygintdiff = 0;
	ybintdiff = 0;
	ig = 1;
	ib = 1;
	i = 1;
	gtrust = 1.0;
	btrust = 0.0;
	baroready = false;
	gnssready = false;
	yb = yybaro(ib);
	yg = yygnss(ig);

	yygtrust = [];
	yybtrust = [];
	yygdiffint = [];
	yybdiffint = [];
	ttbarocorr = [];
	yybarocorr = [];
	ttlag = [];
	yygnsslag = [];
	yybarolag = [];
	% Fused height
	ttfuse = [];
	yyfuse = [];

	while ~isnan(ib) && ~isnan(ig)
		% Previous sensor measurement
		yb = yybaro(ib);
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
			dt = ttgnss(ig) - ttgnss(ig - 1);
			ydiff = abs(yg - yygnss(ig));
			yg = yygnss(ig);
			t = ttgnss(ig);
		else
			baroready = true;
			dt = ttbaro(ib) - ttbaro(ib - 1);
			ydiff = abs(yb - yybaro(ib));
			yb = yybaro(ib);
			t = ttbaro(ib);
		end

		% Save corrected barometer values
		yb = yb + barooffset;
		ttbarocorr = [ttbarocorr, t];
		yybarocorr = [yybarocorr, yb];
		if usegnss
			% Update baro offset, apply correction with each GNSS sensor update
			%barooffsetdelta = (yg - yb) * gtrust * bcorr_slowdown;
			barooffsetdelta = (yg - yb) * bcorr_slowdown;
			barooffset = barooffset + barooffsetdelta;
		end

		if baroready && gnssready
			% Synchronized. Rebalance scores when got both baro., and GNSS
			% Update lags
			ttlag = lag_update(ttlag, t, windowsize);
			yygnsslag = lag_update(yygnsslag, yg, windowsize);
			yybarolag = lag_update(yybarolag, yb, windowsize);
			baroready = false;
			gnssready = false;
			if numel(ttlag) >= 2
				dt = ttlag(2) - ttlag(1);
				% Update baroint
				babsdiff = abs(diff(yybarolag));
				ybintdiff = clamp(0.0001, inf, ybintdiff - intcooldown * dt + (sum(babsdiff) - bsigma * numel(babsdiff)) * dt);
				% Update gnssint
				gabsdiff = abs(diff(yygnsslag));
				ygintdiff = clamp(0.0001, inf, ygintdiff - intcooldown * dt + (sum(gabsdiff) - gsigma * numel(gabsdiff)) * dt);

				gtrust = clamp(0, 1, ybintdiff / (ybintdiff + ygintdiff));
				btrust = 1 - gtrust;
			end

			% Save fused
			ttfuse = [ttfuse, t];
			yyfuse = [yyfuse, gtrust * yg + btrust * yb];
			% Save trust scores
			yybtrust = [yybtrust, btrust];
			yygtrust = [yygtrust, gtrust];
			% Save diffints
			yybdiffint = [yybdiffint, ybintdiff];
			yygdiffint = [yygdiffint, ygintdiff];
		end

		i = i + 1
	end

	save cool;
	figure;
	hold on
	plot(ttbaro, yybaro + ttyygnss(1, 5) - yybaro(1));
	plot(ttgnss, yygnss);
	plot(ttbarocorr, yybarocorr);
	plot(ttfuse, yyfuse);
	legend(
		'baro',
		'gnss',
		'baro corrected',
		'fused'
	);

	% Plot trust scores
	figure;
	hold on;
	plot(ttfuse, yybtrust);
	plot(ttfuse, yygtrust);
	legend("Baro trust", "GNSS trust");

	% Plot integrated diffs
	figure;
	hold on;
	plot(ttfuse, yybdiffint);
	plot(ttfuse, yygdiffint);
	legend("Baro diff int", "GNSS diff int")

	% Plot alongside
	figure;
	ax1 = subplot(2, 1, 1);
	hold on
	plot(ttbaro, yybaro + ttyygnss(1, 5) - yybaro(1));
	plot(ttgnss, yygnss);
	plot(ttbarocorr, yybarocorr);
	plot(ttfuse, yyfuse);
	legend('baro', 'gnss', 'baro corrected', 'fused');
	% Plot integrated diffs
	ax2 =subplot(2, 1, 2);
	hold on;
	plot(ttfuse, yybdiffint);
	plot(ttfuse, yygdiffint);
	legend("Baro diff int", "GNSS diff int")
	linkaxes([ax1, ax2], "x")
end

function [v] = clamp(a, b, v)
	assert(a < b);
	v = min(b, max(a, v));
end

function [lag] = lag_update(lag, val, winsize)
	lag = [lag, val];
	if numel(lag) > winsize
		lag = lag(2:end);
		lag = lag(1:winsize);
	end
end

function [slice] = get_backlag(i, windowsize, series)
	i0 = max([1, i - windowsize]);
	slice = series(i0:i);
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
