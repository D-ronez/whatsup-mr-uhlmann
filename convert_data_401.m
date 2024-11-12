%
% Convert data from raw sensors on GS 401 into actual values
%

function [ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401(path, sid, gnssnoisevar)
	% Extract data from GS log files named like "session_00284_IntStaticPressure.csv"
	% path: path prefix
	% sid: session id
	% ttyygnss: matrix [TIME, LON, LAT, ALT] (s, deg, deg, m)
	% ttybaroalt: matrix [TIME, ALTITUDE] (s m)
	% ttyyaccgyro: matrix [TIME ACCX ACCY ACCZ TEMP GYROX GYROY GYROZ] (s m/s2 m/s2 m/s2 ~ deg/s deg/s deg/s)
	% gnssnoisevar: Additional noise for GNSS data, because is rounded to 1m, not usable for real applications

	pressname = [path "session_" sid "_" "Altitude.csv"];
	ttyybaroalt = csvread(pressname)

	gnssname = [path "session_" sid "_" "NavPosition.csv"];
	ttygnss = csvread(gnssname);

	accgyroname = [path "session_" sid "_" "RawAccelGyroData.csv"];
	ttyyaccgyroraw = csvread(accgyroname);
	ttyaccgyro = process_mpu60xx_accgyro_data(ttyaccgyroraw)
end

function [ttyybaroalt] = process_pressure_data(ttyypressure)
	% ttyypressure: matrix [TIME, PRESSURE] (s Pa)
	% ttybaroalt: matrix [TIME, ALTITUDE] (s m)
	% Converts barometric pressure into barometric altitude. Uses first point as zero reference.
	groundpressure = mean(ttyypressure(1:3, 2)) % Just a mean from first values
	ttybaroalt = ttypressure
	for i in size(ttypressure(1))
			val = ttypressure(i, 2)
			tmp = (val / groundpressure) ^ (1 / 5.255)
			ttybaroalt = 44330 * (1 - tmp)
	end
end

function [ttyyaccgyro] = process_mpu60xx_accgyro_data(ttyaccgyro_raw)
	% Converts raw sensor values from PMU60xx into accelerometer, and gyroscope values
	% Acc/Gyro scaling values
	% ttyaccgyro_raw: matrix [TIME ACCXRAW ACCYRAW ACCZRAW TEMP GYROXRAW GYROYRAW GYROZRAW] (s ADC ADC ADC ~ ADC ADC ADC)
	% ttyyaccgyro: matrix [TIME ACCX ACCY ACCZ TEMP GYROX GYROY GYROZ] (s m/s2 m/s2 m/s2 ~ deg/s deg/s deg/s)
	g = 9.80665
	accscale = 16 % [g]
	gyroscale = 2000 % [deg/s]
	resolution = 2^16 - 1 % 16 bit ADC

	ttyyaccgyro = ttyaccgyro_raw
	ttyaccgyro = ttyaccgyro_raw(:, 2) / resolution * accscale * g
	ttyaccgyro = ttyaccgyro_raw(:, 3) / resolution * accscale * g
	ttyaccgyro = ttyaccgyro_raw(:, 4) / resolution * accscale * g
	ttyaccgyro = ttyaccgyro_raw(:, 6) / resolution * gyroscale
	ttyaccgyro = ttyaccgyro_raw(:, 7) / resolution * gyroscale
	ttyaccgyro = ttyaccgyro_raw(:, 8) / resolution * gyroscale
end
