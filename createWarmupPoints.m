function warmUp = createWarmupPoints(model,bounds,npts)
if nargin<3 || npts<2*size(model.S,1)
    %default # of warmup points
    npts = 2*size(model.S,1);
end

%set concentration bounds
if nargin<2
    %problem setup with slack variables for all inequalities
    bounds = setupMetLP(model);     
else
    
    %setup the initial LP problem to obtain cocnetrations
    
%     if isfield(bounds,'lb')
%         lb = bounds.lb;
%         bndflag_lb = 0;
%     end
%     if isfield(bounds,'ub')
%         ub = bounds.ub;
%         bndflag_ub = 0;
%     end
end

if npts<2*size(bounds.S,1);
     %change deafult to # sampled metabolites only   
    npts = 2*size(bounds.S,1);
end

%not all metabolites are sampled - rearranging lb and ub for all
%metabollites
lb = separate_slack(bounds.lb,model,bounds);
ub = separate_slack(bounds.ub,model,bounds);

%total sampled metabolites
n_mets = length(bounds.mets);
% n_mets = size(bounds.A,2);

%total metabolites
% nmets = size(model.S,1);
warmUp = zeros(n_mets,npts);

ipt = 1;
%generate points
while ipt<=npts/2
    validflag = 0;
%     fprintf('Max/Min metabolite %d %s\n',ipt,bounds.mets{ipt});
    %create random objective function    
%     prxnid = 0;
    %objective for sampled mets only
    bounds.cprod = rand(1,n_mets)-0.5;
    
    %max/min objective
    if ipt<=n_mets
%         prxnid = ipt;
        %objective for sampled mets only
        bounds.cprod = sparse(1,ipt,1,1,n_mets);                
    end
    
    %get maximum and minimum for cprod(ipt)
    [LPmax,LPmin] = solvemetLP(bounds);
    
    %use max points
    if LPmax.flag>0
        xmax = separate_slack(LPmax.x,model,bounds);
        validflag = validflag+1;
    else
        validflag = validflag-1;
    end
    if LPmin.flag>0
        xmin = separate_slack(LPmin.x,model,bounds);
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
            fprintf('%d. %s',ipt,bounds.mets{ipt});
        elseif validflag == 0
            fprintf('Either min or max failed\n');
            fprintf('%d. %s MaxFlag %d MinFlag %d\n',...
                    ipt,bounds.mets{ipt},LPmax.flag,LPmin.flag);
        end
    end  
end

centrepoint = mean(warmUp,2);

%moving in points
warmUp = warmUp*.99+.01*centrepoint*ones(1,npts);

%lb and ub are for lnP and lnS
% warmUp = exp(warmUp);

% warmUp = assignConc(warmUp,model,bounds);
            