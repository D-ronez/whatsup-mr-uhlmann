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
