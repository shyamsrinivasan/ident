function createWarmupPoints(model,bounds,npts)
bndflag_lb = 1;
bndflag_ub = 1;
if nargin<3
    npts = 500;
end
if nargin<2
    lb = model.lb;
    ub = model.ub;
else
    if isfield(bounds,'lb')
        lb = bounds.lb;
        bndflag_lb = 0;
    end
    if isfield(bounds,'ub')
        ub = bounds.ub;
        bndflag_ub = 0;
    end
end

%default # of warmup points
[nmets,nrxns] = size(model.S);
if npts<2*nmets
    npts = 2*nmets;
end

%set reaction bounds if not set
if ~bndflag_ub || ~bndflag_lb
    bounds = setMETbounds(model);
    lb = bounds.lb;
    ub = bounds.ub;
end

warmUp = sparse(nmets,npts);
model = setupMetLP(model);

validflag = 0;
ipt = 1;
%generate points
while ipt<npts/2
    %create random objective function
    bounds.cprod = rand(1,nmets)-0.5;
    
    %max/min objective
    if ipt<=nmets
        bounds.cprod = sparse(1,ipt,1,1,nmets);
    end
    
    %get maximum and minimum for cprod(ipt)
    [LPmax,LPmin] = solvemetLP(model,bounds);
    
    %use max points
    if LPmax.flag>0
        xmax = LPmax.x;
        validflag = validflag+1;
    else
        validflag = validflag-1;
    end
    if LPmin.flag>0
        xmin = LPmin.x;
        validflag = validflag+1;
    else
        validflag = validflag-1;
    end
    
    %move points to within bounds
    xmax(xmax>ub) = ub(xmax>ub);
    xmax(xmax<lb) = lb(xmax<lb);
    
    xmin(xmin>ub) = ub(xmin>ub);
    xmin(xmin<lb) = lb(xmin<lb);
    
    %store points
    warmUp(:,2*ipt-1) = xmin;
    warmUp(:,2*ipt) = xmax;
    
    if validflag == 2
        ipt = ipt+1;
    else
        if validflag<0
            fprintf('Both min and max failed\n');
        elseif validflag == 0
            fprintf('Either min or max failed\n');
        end
    end  
end
            