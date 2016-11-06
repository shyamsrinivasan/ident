function fixed_points = kotte_branches(npoints,range,contdir,eqpts,model,pvec,ap,opts)

% set parameters
global c
c = 3;

neqpts = size(eqpts,1);

% initialize parameter
delpar = (range(2) - range(1))/(npoints-1);

% direction of continuation
if contdir == 1
    pararray = range(1):delpar:range(2);
elseif contdir == -1
    pararray = range(2):-delpar:range(1);
end

arraycont = zeros(neqpts,3*npoints);
% arraycont(:,1:3) = arraycont;

conttype = zeros(neqpts,npoints);

% continuation of fixed points
contpt = continuation(model,pvec,pararray,ap,eqpts,opts);
for i = 1:npoints
%     fixed_points = 
end