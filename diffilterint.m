% Integral differential filter
%

function [x, istate] = difffilter(z, i, istate)
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

	barooffset = 0; % The baro values will be constantly corrected against GNSS
	bcorr_slowdown = 0.008; % Correction scale for baro offset, the lower it is, the slower baro values are corrected
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3) + barooffset;;
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);

	% Fused height
	ttfuse = [];
	yyfuse = [];
	% Corrected baro values
	ttbarocorr = [];
	yybarocorr = [];

	ig = 1;
	ib = 1;
	i = 1;
	gtrust = 1.0;
	btrust = 0.0;
	yb = yybaro(ib);
	yg = yygnss(ig);
	% Difference values used for calculating the score
	trustlag = 0.1;

	yygtrust = [];
	yybtrust = [];

	% DiffLag-s are queues of absolute differences between consecutive sensor measurements
	bdifflag = [];
	gdifflag = [];
	gwinsize = 5;
	bwinsize = gwinsize;

	% For time synchronization. Pre: sensor values' streams are continuous, and do not abrupt
	gnssready = true;
	baroready = true;
	should_rebalance = true;

	btrustint = 0;
	gtrustint = 0;

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
			ygdiff = abs(yg - yygnss(ig));
			% Time synchronization
			if gnssready
				should_rebalance = true;
				gnssready = false;
				baroready = true;
				gdifflag = lag_update(gdifflag, ygdiff, gwinsize);
			end
			yg = yygnss(ig);
			t = ttgnss(ig);
		else
			ybdiff = abs(yb - yybaro(ib));
			% Time synchronization
			if baroready
				baroready = false;
				gnssready = true;
				bdifflag = lag_update(bdifflag, ybdiff, bwinsize);
				should_rebalance = true;
			end
			yb = yybaro(ib);
			t = ttbaro(ib);
		end

		% Rebalance scores
		if should_rebalance
			should_rebalance = false
			newgtrust = clamp(0, 1, sum(bdifflag) / (sum(bdifflag) + sum(gdifflag) + 0.0001));
			newbtrust = 1 - newgtrust;

			gtrustint = gtrustint + sum(bdifflag) * gtrustint;
			btrustint = btrustint + sum(gdifflag) * btrustint;

			gtrust = clamp(0.0001, 1, gtrustint / (gtrustint + btrustint));
			btrust = 1 - gtrust;
		end

		% Save corrected barometer values
		yb = yb + barooffset;
		ttbarocorr = [ttbarocorr, t];
		yybarocorr = [yybarocorr, yb];
		if usegnss
			% Update baro offset, apply correction with each GNSS sensor update
			barooffsetdelta = (yg - yb) * gtrust * bcorr_slowdown;
			barooffset = barooffset + barooffsetdelta;
		end

		% Save fused
		ttfuse(i) = t;
		yyfuse(i) = gtrust * yg + btrust * yb;

		% Save trust scores
		yybtrust = [yybtrust, btrust];
		yygtrust = [yygtrust, gtrust];

		i = i + 1
	end

	save varf;
	figure;
	hold on
	plot(ttbaro, yybaro + ttyygnss(1, 5) - yybaro(1));
	plot(ttgnss, yygnss, 'LineWidth', 4);
	plot(ttbarocorr, yybarocorr);
	plot(ttfuse, yyfuse, 'LineWidth', 3);
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