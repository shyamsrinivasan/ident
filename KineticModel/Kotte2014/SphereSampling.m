% sample npts points from a ball of radius Rs around any point p
function sampledpts = SphereSampling(origin,Rs,npts)
if nargin<3
    npts = 100;
end
if nargin<2
    % ball of unit radius in 3d cartesian space
    Rs = 1;
end
if nargin<1
    origin = [0;0;0];
end

nvar = 3;
% get npts uniformly distributed points
U = random('Uniform',0,1,[nvar npts]);
% get nptsx3 standard normal points
x = random('Normal',0,1,[nvar npts]);
% calculate direction vector
d = x./repmat(sqrt(x(1,:).^2+x(2,:).^2+x(3,:).^2),nvar,1);
% calculate distance from origin
lambda = Rs*U.^(1/3);

% translation points to new origin in origin
sampledpts = repmat(origin,1,npts)+lambda.*d;