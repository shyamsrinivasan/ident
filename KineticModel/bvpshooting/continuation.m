% continuation for 2 point BVP
% define problem
% Holt.m
% initial conditions
nvar = 5;
yinit = zeros(nvar,1);

% known initial conditios @ t0
yiknwn = zeros(nvar,1);
yiknwn(1) = 1; yiknwn(2) = 1; yiknwn(4) = 1;
% terminal conditions
yterm = zeros(nvar,1);
yterm(4) = 1;

yiknwn = find(yiknwn);
r = length(yiknwn);
yiunkwn = setdiff(1:nvar,yiknwn);

% known terminal conditions @ tf
yfknwn = zeros(nvar,1);
yfknwn(2) = 1; yfknwn(4) = 1;
yfunkwn = find(~yfknwn);
yfknwn = find(yfknwn);

% choose unknown initial conditions
yinit(yiunkwn) = [-1;0.6];

% integrate till time t1 till which there are no numerical problems
% solve 2 point bvp over (t0,t1) using goodman lance method
% extend time interval from (t0,t1) to (t0,t2)