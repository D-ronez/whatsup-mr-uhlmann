42;

function [rot] = imu_rot_update(rot, a, g, dt)
	% Updates UAV's rotation matrix using new IMU values. Implements ~~complementary~~ filter
end

function [rot] = imu_rotation(gvec)
	% Estimates the rotation matrix from accelerometer readings. Assumes that the only force acting on the UAV
	% is gravitational force (holds true, when the vehicle is stationary)
	% gvec - measurement of the accelerometer
	rot = zeros(3, 3);
	rot(:, 3) = normalize(gvec);
	% Given the absence of other data (magnetometer), UAV's X is residing in global YZ plane
	rot(:, 1) = normalize(crossProduct3d(rot(:, 3)', [0 1 0]));
	rot(:, 2) = normalize(crossProduct3d(rot(:, 3)', rot(:, 1)'));
end
