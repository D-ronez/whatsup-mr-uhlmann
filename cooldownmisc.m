% Imitates continuous flow of data

function [] = cooldownmisc(filtertype)
	[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

	barooffset = ttyygnss(1, 1); % The baro values will be constantly corrected against GNSS
	bcorr_slowdown = 0.008; % Correction scale for baro offset, the lower it is, the slower baro values are corrected
	ttbaro = ttyybaroalt(:, 1)';
	yybaro = ttyybaroalt(:, 3)' + barooffset;
	ttgnss = ttyygnss(:, 1)';
	yygnss = ttyygnss(:, 5)';

	% Cooldown for accumulated error
	msewinsize = 100;
	maewinsize = 100;

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
	yygint = [];
	yybint = [];
	ttbarocorr = [];
	yybarocorr = [];
	% Fused height
	ttfuse = [];
	yyfuse = [];
	% MSE
	ttmse = [];
	yymse = [];
	% MAE
	ttmae = [];
	yymae = [];
	% GNSS for MSE, MAE estimations to preserve dimensionality w/ fused values
	yygnssmxe = [];

	on_baro = nan;
	on_gnss = nan;
	if strcmp(filtertype, 'var')
		on_baro = @(bt, gt, t, y)varcool_on_baro(bt, gt, t ,y);
		on_gnss = @(bt, gt, t, y)varcool_on_gnss(bt, gt, t ,y);
		init = @()varcool_init();
	else
		error('wrong filter type');
	end

	init();
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
			t = ttgnss(ig);
			gnssready = true;
		else
			yb = yybaro(ib);
			t = ttbaro(ib);
			baroready = true;
		end

		% Save corrected barometer values
		yb = yb + barooffset;
		ttbarocorr = [ttbarocorr, t];
		yybarocorr = [yybarocorr, yb];
		if usegnss
			% Update baro offset, apply correction with each GNSS sensor update
			barooffsetdelta = (yg - yb) * bcorr_slowdown;
			barooffset = barooffset + barooffsetdelta;
		end

		% Update trust scores
		if usegnss
			[btrust, gtrust] = on_gnss(btrust, gtrust, t, yg);
		else
			[btrust, gtrust] = on_baro(btrust, gtrust, t, yb);
		end

		% New estimation, when both sensors are ready
		if gnssready && baroready
			global ybint;
			global yybint;
			global ygint;
			global yygint;
			baroready = false;
			gnssready = false;

			% Save fused
			ttfuse = [ttfuse, t];
			yyfuse = [yyfuse, gtrust * yg + btrust * yb];
			% Save GNSS for MAE, MSE
			yygnssmxe = [yygnssmxe, yg];
			% Save trust scores
			yybtrust = [yybtrust, btrust];
			yygtrust = [yygtrust, gtrust];
			% Save ints
			yybint = [yybint, ybint];
			yygint = [yygint, ygint];
			% Save MSE
			winend = numel(yyfuse);
			winbegin = max(1, winend - msewinsize);
			a = yyfuse(winbegin:winend);
			b = yygnssmxe(winbegin:winend);
			ttmse = [ttmse, t];
			yymse = [yymse, sum(a - b) .^ 2 / msewinsize];
			% Save MAE
			winend = numel(yyfuse);
			winbegin = max(1, winend - maewinsize);
			a = yyfuse(winbegin:winend);
			b = yygnssmxe(winbegin:winend);
			ttmae = [ttmae, t];
			yymae = [yymae, sum(abs(a - b)) / maewinsize];
		end

		i = i + 1
	end

	save coolmisc;
	figure;
	hold on
	plot(ttbaro, yybaro + ttyygnss(1, 5) - yybaro(1));
	plot(ttgnss, yygnss);
	plot(ttfuse, yyfuse);
	legend(
		'baro',
		'gnss',
		'fused'
	);

	% Plot trust scores
	if false
		figure;
		hold on;
		plot(ttfuse, yybtrust);
		plot(ttfuse, yygtrust);
		legend("Baro trust", "GNSS trust");
	end

	% Plot integrated diffs
	% figure;
	% hold on;
	% plot(ttfuse, yybint);
	% plot(ttfuse, yygint);
	% legend("Baro diff int", "GNSS diff int")

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
	plot(ttfuse, yybint);
	plot(ttfuse, yygint);
	legend("Baro diff int", "GNSS diff int")
	linkaxes([ax1, ax2], "x")

	% Plot MSE
	figure;
	plot(ttmse, yymse);

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

source("varslide.m")

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

function [] = varcool_init()
	source("coolglobal.m")
	yybarolag = [];
	ttbarolag = [];
	ttgnsslag = [];
	yygnsslag = [];
	yygint = [];
	yybint = [];
	baroready = false;
	gnssready = false;
	ybint = 0;
	ygint = 0;
	tprev = 0;
	intcooldown = 0.2; % [m/s]
	bsigma = 0.5; % [m]
	gsigma = 0.1; % [m]
	windowsize = 15;
end

function [btrust, gtrust] = varcool_on_baro(btrust, gtrust, t, y)
	source("coolglobal.m")
	baroready = true;
	yb = y;
	[btrust, gtrust] = varcool(btrust, gtrust, t);
end

function [btrust, gtrust] = varcool_on_gnss(btrust, gtrust, t, y)
	source("coolglobal.m")
	gnssready = true;
	yg = y;
	[btrust, gtrust] = varcool(btrust, gtrust, t);
end

function [btrust, gtrust] = varcool(btrust, gtrust, t)
	source("coolglobal.m")
	if baroready && gnssready
		baroready = false;
		gnssready = false;
		if isnan(tprev)
			tprev = t;
			return;
		end
		dt = t - tprev;
		% Synchronized. Rebalance scores when got both baro., and GNSS
		% Update lags
		yygnsslag = lag_update(yygnsslag, yg, windowsize);
		yybarolag = lag_update(yybarolag, yb, windowsize);
		baroready = false;
		gnssready = false;

		backlag = yybarolag;
		ybint = clamp(0.0001, inf, ybint - intcooldown * dt + (std(backlag, 1) - bsigma) * dt);
		backlag = yygnsslag;
		ygint = clamp(0.0001, inf, ygint - intcooldown * dt + (std(backlag, 1) - gsigma) * dt);
		gtrust = clamp(0, 1, ybint / (ybint + ygint));
		btrust = 1 - gtrust;
		tprev = t;
	end
end

