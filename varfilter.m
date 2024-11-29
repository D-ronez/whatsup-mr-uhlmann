% Acquires fused altitude value from baro, and GNSS, and performs scaled
% on-the-fly correction of baro using GNSS to prevent it from drifting away over
% time.
%
% Variance-based filter for generating fused altitude from GNSS, and baro
% measurements. It computes trust scores for GNSS, and baro, and produces fused
% value based on weighted sum of both sensors whose weights are inferred
% according to their respective trust values.
%

function [x, istate] = varfilter(z, i, istate)
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

	barooffset = 0; % The baro values will be constantly corrected against GNSS
	bcorr_slowdown = 0.008; % Correction scale for baro offset, the lower it is, the slower baro values are corrected
	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3) + barooffset;;
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);

	% Lag values are used for calculating trust scores
	glag = [yygnss(1)];
	blag = [ttgnss(1)];
	gwinsize = 10; % Size of history queue for GNSS values
	bwinsize = gwinsize * 5;

	% Fused height
	ttfuse = [];
	yyfuse = [];
	% Corrected baro values
	ttbarocorr = [];
	yybarocorr = [];

	ig = 1;
	ib = 1;
	i = 1;
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
			yg = yygnss(ig);
			glag = lag_update(glag, yg, gwinsize);
			t = ttgnss(ig);
		else
			yb = yybaro(ib);
			blag = lag_update(blag, yb, bwinsize);
			t = ttbaro(ib);
		end

		% Save corrected barometer values
		yb = yb + barooffset;
		ttbarocorr = [ttbarocorr, t];
		yybarocorr = [yybarocorr, yb];

		% Calculate trust scores
		vb = var(blag);
		vg = var(glag);
		gtrust = vb / (vb + vg);
		btrust = 1 - gtrust;

		if usegnss
			% Update baro offset, apply correction with each GNSS sensor update
			barooffsetdelta = (yg - yb) * gtrust * bcorr_slowdown;
			barooffset = barooffset + barooffsetdelta;
		end

		% Save fused
		ttfuse(i) = t;
		yyfuse(i) = gtrust * yg + btrust * yb;

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
end

function [v] = clamp(a, b, v)
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
