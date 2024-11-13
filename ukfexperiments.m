42; % This is script

pkg load matgeom;
source("util.m")
source("imu_estimate.m")

[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");

% Imu calculations
% Z vector in global frame. The UAV is stable, acc. measurements are aligned with the gravitational vector (G-vector)
gvec = [ttyyaccgyro(1, 3), ttyyaccgyro(1, 4), ttyyaccgyro(1, 5)];
rot = imu_rotation_from_acc(gvec);
figure
grid on
hold on
plot3_rot(rot)
[roll, pitch, yaw] = rot2euler(rot);
yaw = 0
hold on;
plot3_euler([roll pitch yaw]);

% Estimate orientation changes

if false %
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
end %
