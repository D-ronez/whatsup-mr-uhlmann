%
% Convert data from raw sensors on GS 401 into actual values
%

function [ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401(path, sid, gnssnoisevar)
	% Extract data from GS log files named like "session_00284_IntStaticPressure.csv"
	% path: path prefix
	% sid: session id
	% ttyygnss: matrix [TIME, NAME, LON, LAT, ALT] (s, deg, deg, m)
	% ttybaroalt: matrix [TIME, NAME, ALTITUDE] (s m)
	% ttyyaccgyro: matrix [TIME NAME ACCX ACCY ACCZ TEMP GYROX GYROY GYROZ] (s m/s2 m/s2 m/s2 ~ deg/s deg/s deg/s)
	% gnssnoisevar: Additional noise for GNSS data, because is rounded to 1m, not usable for real applications

	pressname = [path "session_" sid "_" "Altitude.csv"];
	ttyybaroalt = csvread(pressname);

	gnssname = [path "session_" sid "_" "NavPosition.csv"];
	ttyygnss = csvread(gnssname);

	accgyroname = [path "session_" sid "_" "RawAccelGyroData.csv"];
	ttyyaccgyroraw = csvread(accgyroname);
	ttyyaccgyro = process_mpu60xx_accgyro_data(ttyyaccgyroraw);
end

function [ttyybaroalt] = process_pressure_data(ttyypressure)
	% ttyypressure: matrix [TIME, NAME, PRESSURE] (s Pa)
	% ttybaroalt: matrix [TIME, NAME, ALTITUDE] (s m)
	% Converts barometric pressure into barometric altitude. Uses first point as zero reference.
	groundpressure = mean(ttyypressure(1:3, 2)) % Just a mean from first values
	ttyybaroalt = ttyypressure
	for i = 1:size(ttyypressure(1))
			val = ttyypressure(i, 2);
			tmp = (val / groundpressure) ^ (1 / 5.255);
			ttybaroalt = 44330 * (1 - tmp);
	end
end

function [ttyyaccgyro] = process_mpu60xx_accgyro_data(ttyaccgyro_raw)
	% Converts raw sensor values from PMU60xx into accelerometer, and gyroscope values
	% Acc/Gyro scaling values
	% ttyaccgyro_raw: matrix [TIME NAME ACCXRAW ACCYRAW ACCZRAW TEMP GYROXRAW GYROYRAW GYROZRAW] (s ADC ADC ADC ~ ADC ADC ADC)
	% ttyyaccgyro: matrix [TIME NAME ACCX ACCY ACCZ TEMP GYROX GYROY GYROZ] (s m/s2 m/s2 m/s2 ~ deg/s deg/s deg/s)
	g = 9.80665;
	accscale = 16; % [g]
	gyroscale = 2000; % [deg/s]
	resolution = 2^16 - 1; % 16 bit ADC

	ttyyaccgyro = ttyaccgyro_raw;
	ttyyaccgyro(:, 3) = ttyaccgyro_raw(:, 3) / resolution * accscale * g;
	ttyyaccgyro(:, 4) = ttyaccgyro_raw(:, 4) / resolution * accscale * g;
	ttyyaccgyro(:, 5) = ttyaccgyro_raw(:, 5) / resolution * accscale * g;
	ttyyaccgyro(:, 7) = ttyaccgyro_raw(:, 7) / resolution * gyroscale;
	ttyyaccgyro(:, 8) = ttyaccgyro_raw(:, 8) / resolution * gyroscale;
	ttyyaccgyro(:, 9) = ttyaccgyro_raw(:, 9) / resolution * gyroscale;
end
