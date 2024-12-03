if true
	load ukfgnss;
	ttmseukf = ttmse;
	yymseukf = yymse;
	ttfuseukf = ttfuse;
	yyfuseukf = yyfuse;
	ttmaeukf = ttmae;
	yymaeukf = yymae;
	load coolvar;
	ttmsecoolv = ttmse;
	yymsecoolv = yymse;
	ttfusecoolv = ttfuse;
	yyfusecoolv = yyfuse;
	ttmaecoolv = ttmae;
	yymaecoolv = yymae;

	% MSE
	figure
	hold on
	title("MSE with sliding window")
	plot(ttmseukf, yymseukf);
	plot(ttmsecoolv, yymsecoolv);
	legend("Sliding window MSE UKF", "Sliding window MSE Var-Cooldown Filter");

	% MAE
	figure
	hold on
	title("MAE with sliding window")
	plot(ttmaeukf, yymaeukf);
	plot(ttmaecoolv, yymaecoolv);
	legend("Sliding window MAE UKF", "Sliding window mae Var-Cooldown Filter");

	% Results comparison
	figure;
	hold on;
	plot(ttgnss, yygnss, 'Color', 'black');
	plot(ttfusecoolv, yyfusecoolv, 'Color', 'red');
	plot(ttfuseukf, yyfuseukf, 'Color', 'blue');
	legend("GNSS", "Cooldown filter", "UKF");

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
