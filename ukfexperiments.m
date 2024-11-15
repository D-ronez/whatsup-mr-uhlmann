42; % This is script

pkg load matgeom;
source("util.m")
source("imu_estimate.m")
addpath("ukf")

[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

% Imu calculations
% Z vector in global frame. The UAV is stable, acc. measurements are aligned with the gravitational vector (G-vector)
gvec = [ttyyaccgyro(1, 3), ttyyaccgyro(1, 4), ttyyaccgyro(1, 5)];
rot = imu_rotation_from_acc(gvec);
grid on
hold on
%plot3_rot(rot)
euler = rot2euler(rot);
yaw = 0
hold on;

if false % Plot vertial acceleration graph along w/ baro alt
	tt = zeros(size(ttyyaccgyro)(1), 1);
	yy = zeros(size(ttyyaccgyro)(1), 1);
	t = ttyyaccgyro(1, 1);
	tt(1) = t;
	for i = 2:size(ttyyaccgyro)(1)
		% Prepare data
		acc = ttyyaccgyro(i, 3:5);
		gyro = ttyyaccgyro(i, 7:9);
		prevt = t;
		t = ttyyaccgyro(i, 1);
		dt = t - prevt;
		% Update rotation matrix in global frame
		euler = euler_estimate_complemetary(euler, acc, gyro, dt);
		i / size(ttyyaccgyro)(1)
		% Get gravitational component
		az = vert_accel(euler, acc');
		% Fix data
		tt(i) = t;
		yy(i) = az;
	end

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
end

if false % Animated plot of position estimation
	startpos = 1
	printeach = 10
	updateeach = printeach * 5
	t = ttyyaccgyro(startpos, 1)
	for i = startpos:size(ttyyaccgyro)(1)
		acc = ttyyaccgyro(i, 3:5);
		gyro = ttyyaccgyro(i, 7:9);
		prevt = t;
		t = ttyyaccgyro(i, 1);
		dt = t - prevt;
		euler = euler_estimate_complemetary(euler, acc, gyro, dt);
		if mod(i, updateeach) == 0
			cla
		end
		if mod(i, printeach) == 0
			hold on;
			xlim([-0.3 0.3])
			ylim([-0.3 0.3])
			zlim([-0.3 0.3])
			plot3_euler(euler);
			pause(0.02)
			t
		end
	end
end

if false % Plot raw data
	hold on;
	grid on;
	% Baro
	ax1 = subplot(2, 1, 1);
	ttbaro = ttyybaroalt(:, 1)
	yybaro = ttyybaroalt(:, 3)
	plot(ttbaro, yybaro);
	% GNSS
	ax2 = subplot(2, 1, 2);
	ttgnss = ttyygnss(:, 1);
	yygnss = ttyygnss(:, 5);
	plot(ttgnss, yygnss)
	linkaxes([ax1, ax2], "x")
	% Acc / Gyro
	figure
	subplot(2, 1, 1)
	hold on
	grid on
	ttacc = ttyyaccgyro(:, 1);
	for i = [3 4 5]
		plot(ttacc, ttyyaccgyro(:, i));
	end
	subplot(2, 1, 2)
	for i = [7 8 9]
		hold on
		grid on
		plot(ttacc, ttyyaccgyro(:, i));
	end
end
