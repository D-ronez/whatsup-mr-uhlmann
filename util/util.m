% Conventions
%
% ROTATION MATRIX (ROT):
% Represents rotation of a vector in another frame of reference:
% rot: rot = [Xdir, Ydir, Zdir]. Xdir - column, direction of x in another frame of
% reference

42;

% For 3d vector plotting
global QUIVER_SCALE = 0.1

function [] = plot3_vector2(pa, pb)
	global QUIVER_SCALE
	pa = pa * QUIVER_SCALE;
	pb = pb * QUIVER_SCALE;
	quiver3(pa(1), pa(2), pa(3), pb(1), pb(2), pb(3))
end

function [] = plot3_rot(rot)
	% rot: rotation matrix
	plot3_vector(rot(:, 1)')
	hold on
	plot3_vector(rot(:, 2)')
	hold on
	plot3_vector(rot(:, 3)')
end

function [rot] = euler2rot(euler)
	roll = euler(1);
	pitch = euler(2);
	yaw = euler(3);

	rollrot = [
		1 0 0;
		0 cos(roll) -sin(roll);
		0 sin(roll) cos(roll);
	];
	pitchrot = [
		cos(pitch) 0 sin(pitch);
		0 1 0;
		-sin(pitch) 0 cos(pitch);
	];
	yawrot = [
		cos(yaw) -sin(yaw) 0;
		sin(yaw) cos(yaw) 0;
		0 0 1;
	];
	rot = rollrot * pitchrot * yawrot;
end

function [] = plot3_euler(euler)
	% Same as plot3_rot, but uses euler angles
	% roll: [rad]
	% pitch: [rad]

	rot = euler2rot(euler);
	plot3_rot(rot);
end

function [] = plot3_vector(p)
	plot3_vector2([0 0 0], p);
end

function euler = rot2euler(rot)
	% returns [roll pitch yaw]
	yaw = atan2(rot(2, 1), rot(2, 2));
	roll = atan2(rot(3, 2), rot(2, 2));
	pitch = atan2(rot(1, 3), rot (1, 1));
	euler = [roll pitch yaw];
end

function [pos] = ts_find_measurement_before(tts, t)
	% Finds measurement before "t" in timeseries' ("tss")

	if size(tts)(2) < 2
		pos = 1;
	end
	for pos = 2:size(tts)
		if t > tts(pos)
			return
		end
	end
end
