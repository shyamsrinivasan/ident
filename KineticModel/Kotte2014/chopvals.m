function [x,y,z,r] = chopvals(x,y,z,chop_pos)

% terminate all trajectories when they cross zero
allx = [];ally = [];allz = [];
for ip = 1:size(x,1)
    % find first negative
    iidx = find(x(ip,:)<0,1,'first');
    iidy = find(y(ip,:)<0,1,'first');
    iidz = find(z(ip,:)<0,1,'first');
    iid = min([iidx iidy iidz]);
    allx = [allx x(ip,1:iid-1)];
    ally = [ally y(ip,1:iid-1)];
    allz = [allz z(ip,1:iid-1)];
end
x = allx;
y = ally;
z = allz;

% remove points past threshold on 3-D plane
r = find(x>chop_pos(1)|x<0|y>chop_pos(2)|y<0|z>chop_pos(3)|z<0);
x(r) = [];
y(r) = [];
z(r) = [];