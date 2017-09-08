function [prob_struct,data] = identopt_setup_v2(data,fixed_pvalue)

if isfield(data,'pname')
    pname = data.pname;
end
% if isfield(data,'casmodelfun')
%     casmodelf = data.casmodelfun;
% end
% if isfield(data,'odefun')
%     intfun = data.integratorfun;
% end
% if isfield(data,'nlaefh')
%     nlaefh = data.nlaefh;
% end
if isfield(data,'objfun')
    objfun = data.objfun;
end
if isfield(data,'xinit')
    xinit = data.xinit;
end
if isfield(data,'xexp')
    xexp = data.xexp;
    yexp = [sum(xexp(1:3,:));sum(xexp(4:6,:));sum(xexp(7:9,:))];
end
if isfield(data,'tspan')
    tspan = data.tspan;
end
% if isfield(data,'freq')
%     freq = data.freq;
% end
% if isfield(data,'p_pert')
%     p_pert = data.p_pert;
% end
% if isfield(data,'p_pert_logical')
%     p_pert_logical = data.p_pert_logical;
% end
npts = length(tspan)-1;
data.npts = npts;

if ~isfield(data,'ident_idx') && ~isempty(pname)    
    plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
            'd','V4max','k1cat','V3max','V2max','acetate'}; 
    ident_idx = find(strcmpi(plist,pname));    
    data.ident_idx = ident_idx;
elseif ~isfield(data,'ident_idx') && isempty(pname)    
    error('No parameter chosen for identifiability analysis');
end    
npar = 13; % length(plist);
np_unchg = 10;

% sort out variables from fixed p values
if ismember(ident_idx,11:npar)
    % constant value could be a perturbation variable
    data.p_pert_logical(end-(npar-ident_idx),:) = 0;    
else
    % constant value is not a perturbation variable
end
nvar = np_unchg+length(find(data.p_pert_logical));
data.np_chang = size(data.p_pert_logical,1);

% objective and bounds
obj = @(xvar)objfun(xvar,fixed_pvalue,data);
[lb,ub] = ident_bounds(nvar);
jacfun = [];
hessfun = [];

% test objfun
scale = ones(10,1);
scale(3) = 1e6;
p_unch = data.odep(1:10)'./scale;
% p_pert = [data.odep(11)';data.odep(13)';data.odep(11)';data.odep(13)'];
p_pert = [2;.1;1;1];
p0 = [p_unch;p_pert];
objval = obj(p0);

prob_struct = struct('objfun',obj,'grad',jacfun,'hess',hessfun,...                
                'lb',lb,'ub',ub,...
                'npts',npts,'xinit',xinit,'xexp',xexp);

