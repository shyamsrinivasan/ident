% planar_circle.m, by Daniel Lyons (March 2014)
% Computes a circle of (npoints+1) points in the x1-x2 plane with radius r.
% Modified only var names and fun name

function [x1, x2] = getplanarcircle(npoints,r)

phase_delta = 2*pi/npoints;
phase = 0:phase_delta:2*pi;

%parametrize unit circle in R^2 (x1,x2)
x1 = r*cos(phase);
x2 = r*sin(phase);
