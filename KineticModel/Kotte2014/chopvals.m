function [x,y,z,r] = chopvals(x,y,z,chop_pos)

% remove points past threshold on 3-D plane
% chop_pos = 5;

r = find(x>chop_pos(1)|x<0|y>chop_pos(2)|y<0|z>chop_pos(3)|z<0);
x(r) = [];
y(r) = [];
z(r) = [];