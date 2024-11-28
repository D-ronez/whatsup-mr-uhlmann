function [x, istate] = varfilter(z, i, istate)
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

	ttbaro = ttyybaroalt(:, 1);
	yybaro = ttyybaroalt(:, 3);
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);

	ig = 1;
	ib = 1;

	glag = [yygnss(1)];
	blag = [ttgnss(1)];

	gwinsize = 10;
	bwinsize = gwinsize * 5;

	ttfuse = [];
	yyfuse = [];
	i = 1;

	while ~isnan(ib) && ~isnan(ig)
		% Previous sensor measurement
		yb = yybaro(ib);
		yg = yygnss(ig);

		[ig, ib, usegnss] = pick_next_point_t(ig, ib, ttgnss, ttbaro);
		t = 0;
		ig
		ib

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

		gtrust = clamp(0, 1, 1 / (var(blag) + 0.0001));
		btrust = 1 - gtrust;

		% Save fused
		ttfuse(i) = t;
		yyfuse(i) = gtrust * yg + btrust * yb;

		i = i + 1
	end

	save varf;
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
	if ia >= numel(tta)
		ia = nan;
	elseif ib >= numel(ttb)
		ib = nan;
	else
		ta = tta(ia + 1);
		tb = ttb(ib + 1);
		if ta > tb
			tb = tb + 1;
			i = 2;
		else
			ta = ta + 1;
			i = 1;
		end
	end
end
