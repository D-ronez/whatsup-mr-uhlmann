1; % This is script

[ttyybaroalt, ttyygnss, ttyyaccgyro] = convert_data_401("/home/dm/Documents/CODE-69910-OverthrottleGnssSag/data/from-pilots/csvs/", "00284");
% Baro data
ttbaro = ttyybaroalt(:, 1)';
yybaro = ttyybaroalt(:, 3)';
% Gnss data
ttgnss = ttyygnss(:, 1);
yygnss = ttyygnss(:, 5);

% Print raw data
figure;
hold on;
grid on;
% Baro
ax1 = subplot(2, 1, 1);
plot(ttbaro, yybaro);
% GNSS
ax2 = subplot(2, 1, 2);
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
