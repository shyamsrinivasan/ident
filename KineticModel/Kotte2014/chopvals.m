function [x,y,z,r] = chopvals(x,y,z,chop_pos)

% remove points past threshold on 3-D plane
% chop_pos = 5;

r = find(x>chop_pos|x<0|y>chop_pos|y<0|z>chop_pos|z<0);
x(r) = [];
y(r) = [];
z(r) = [];