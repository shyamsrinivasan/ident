% planar_circle.m, by Daniel Lyons (March 2014)
% Computes a circle of (npoints+1) points in the x1-x2 plane with radius r.
% Modified only var names and fun name

function [x1,x2,x3] = getplanarcircle(npoints,r)

phase_delta = 2*pi/npoints;
phase = 0:phase_delta:2*pi;
% deltheta = 2*pi/npoints;
% theta = 0:deltheta:2*pi;
% delphi = pi/npoints;
% phi = 0:delphi:pi;
% delu = 2*r/npoints;
% u = -r:delu:r;

% x1 = sqrt(repmat(r^2,1,length(u))-u.^2).*cos(theta);
% x2 = sqrt(repmat(r^2,1,length(u))-u.^2).*sin(theta);
% x3 = u;

%parametrize unit circle in R^2 (x1,x2)
% x1 = r.*cos(theta).*sin(phi);
% x2 = r.*sin(theta).*sin(phi);
% x3 = r.*cos(phi);

x1 = r*cos(phase);
x2 = r*sin(phase);
x3 = [];
