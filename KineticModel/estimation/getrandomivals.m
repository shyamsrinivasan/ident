% choose/supply random initial values for optimization
% eps - ball around given value for sampling
function rnd_pts = getrandomivals(optimdata,eps,nval)
if nargin<3
    nval = 1;
end
if nargin<2
    eps = .05;
end
nvar = optimdata.nvar;
nc = optimdata.nc;
nf = optimdata.nf;
npert = optimdata.npert;

% set initial value
x0 = [optimdata.xexp;optimdata.odep(optimdata.p_id)';optimdata.vexp];

if nval>1
    pd = makedist('Uniform','lower',-eps,'upper',eps);
    sampled_pts = random(pd,nvar,nval);
    rnd_pts = repmat(x0,1,nval).*(1+sampled_pts);
else
    rnd_pts = x0;
end


