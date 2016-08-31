function vectorfield(xdyn,model,pvec)
eps = 1e-4;
sep_start = min(xdyn
npts = size(xdyn,2);
nvar = size(npts,1);
givenModel = @(t,x)KotteODE(t,x,model,pvec);
for ipt = 1:npts    
    zi = saddle+eps*ones(nvar,1);
    zj = saddle-eps*ones(nvar,1);
    Di = givenModel(0,zi);
    Dj = givenModel(0,zj);
end