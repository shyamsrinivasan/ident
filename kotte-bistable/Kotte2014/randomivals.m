function newivals = randomivals(bounds,npts)
nvar = size(bounds,1);
r = zeros(nvar,npts);

% get numbers from a uniform distribution
for i = 1:nvar
    pd = makedist('Uniform','lower',bounds(i,1),'upper',bounds(i,2));
    r1 = random(pd,npts,1);
    r(i,:) = r1';
end

% pd1 = makedist('Uniform','lower',newminx(1),'upper',newmaxx(1));
% pd2 = makedist('Uniform','lower',newminx(2),'upper',newmaxx(2));
% pd3 = makedist('Uniform','lower',newminx(3),'upper',newmaxx(3));
% r1 = random(pd1,ptspinterval,1);
% r2 = random(pd2,ptspinterval,1);
% r3 = random(pd3,ptspinterval,1);

newivals = r';