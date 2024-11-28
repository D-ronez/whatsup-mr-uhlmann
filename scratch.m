
function [] = glb()
	global var = 42
	var
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
