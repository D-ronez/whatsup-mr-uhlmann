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

function [] = plot3_vector(p)
	plot3_vector2([0 0 0], p);
end

function [roll, pitch, yaw] = rot2euler(rot)
	yaw = atan2(rot(2, 1), rot(2, 2));
	roll = atan2(rot(3, 2), rot(2, 2));
	pitch = atan2(rot(1, 3), rot (1, 1));
end
