% Complementary filter, accelerometer coefficient
global COMP_ACC_COEF = 0.9
global COMP_GYRO_COEF = 1.0 - COMP_ACC_COEF

42;

function [rot] = imu_rotation_from_acc(gvec)
	% Estimates the rotation matrix from accelerometer readings. Assumes that the only force acting on the UAV
	% is gravitational force (holds true, when the vehicle is stationary)
	% gvec - measurement of the accelerometer
	% Assumes that the accelerometer is oriented so Z is normal to the ground
	rot = zeros(3, 3);
	rot(:, 3) = normalize(gvec);
	% Given the absence of other data (magnetometer), assume UAV's X is residing
	% in global YZ plane
	rot(:, 1) = normalize(crossProduct3d(rot(:, 3)', [0 1 0]));
	rot(:, 2) = normalize(crossProduct3d(rot(:, 3)', rot(:, 1)'));
end

function [roll, pitch] = imu_euler2_from_acc(gvec)
	% https://seanboe.com/blog/complementary-filters
	gvec = normalize(gvec);
	pitch = atan(-gvec(1) / sqrt(gvec(2)^2 + gvec(3)^2));
	roll = atan(gvec(2) / sqrt(gvec(1)^2 + gvec(3)^2));
end

function [euler] = euler_estimate_complemetary(euler, acc, gyro, dt)
	% Estimates UAV orientation using complementary filter
	% - euler: previous euler estimation. Format [ROLL, PITCH, YAW*]. May have
	% dimensionality 1x2.
	global COMP_ACC_COEF;
	global COMP_GYRO_COEF;

	source("util.m");

	% Assume no forces other than gravity are present
	[aroll, apitch] = imu_euler2_from_acc(acc);
	% Compute Gyro roll, pitch delta estimates
	groll_delta = gyro(1) * dt;
	gpitch_delta = gyro(2) * dt;

	% Calculate the new angle
	roll = euler(1);
	pitch = euler(2);
	roll = aroll * COMP_ACC_COEF + (roll + groll_delta) * COMP_GYRO_COEF;
	pitch = apitch * COMP_ACC_COEF + (pitch + gpitch_delta) * COMP_GYRO_COEF;
	% Update the matrix
	euler(1) = roll;
	euler(2) = pitch;
end
