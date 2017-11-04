% choose/supply random initial values for optimization
% eps - ball around given value for sampling
function rnd_pts = getrandomivals(data,eps,nval)
if nargin<3
    nval = 1;
end
if nargin<2
    eps = .05;
end
if isfield(data,'type')
    type = data.type;
else
    type=1;
end
nvar = data.nvar;
nc = data.nc;
% nf = data.nf;
np = data.np;
npert = data.npert;

% set initial value
if type==1
    x0 = [data.xexp;data.odep(data.p_id)';data.vexp];
elseif type==2
    x0 = [data.xexp;data.odep(data.p_id)';data.vexp;.01;.01];
end
x0(nc*npert+1:nc*npert+np) = .02;

if nval>1
    pd = makedist('Uniform','lower',-eps,'upper',eps);
    sampled_pts = random(pd,nvar,nval);
    rnd_pts = repmat(x0,1,nval).*(1+sampled_pts);
else
    rnd_pts = x0;
end


