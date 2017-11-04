% setup optimization problem for each flux in optimdata input
function [prob,newdata] = setup_prob_typeb(optimdata)

if isfield(optimdata,'flxid')
    flxid = optimdata.flxid;
end
% choose problem type
if isfield(optimdata,'type')
    type = optimdata.type;
end
if isempty(type)
    % type a the original problem where noise is not a var
    type = 1;
end

% all parameters available
rxn_plist = {{'K1ac', 'k1cat'},{},{'K3fdp','K3pep','V3max'},{'K2pep','V2max'},...
              {'V4max'},{}};
relevant_list = rxn_plist{flxid};       
% all parameters available (ordered)
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA','acetate'};        
% find relevant parametyers in ordered list for id
if ~isempty(relevant_list)
    p_id = cellfun(@(x)strcmpi(x,plist),relevant_list,'UniformOutput',false);
    p_id = cellfun(@(x)find(x),p_id);
else    
    error('No parameters to estimate');
end
np = length(p_id);
if isfield(optimdata,'nc')
    nc = optimdata.nc;
end
if isfield(optimdata,'nf')
    nf = optimdata.nf;
end
if isfield(optimdata,'vexp')
    vss_exp = optimdata.vexp(flxid,:);    
end
if isfield(optimdata,'xexp')
    xss_exp = optimdata.xexp;
end
npert = size(xss_exp,2);
vss_exp_v = reshape(vss_exp,[nf*npert,1]);
xss_exp_v = reshape(xss_exp,[nc*npert,1]);

% determine nvar
if type==1
    nvar = nc*npert+nf*npert+np;
elseif type==2
    % 2 noise vars for conc and fluxes
    nvar = nc*npert+nf*npert+np+2;
end
newdata = optimdata;
newdata.vexp = vss_exp_v;
newdata.xexp = xss_exp_v;
newdata.npert = npert;
newdata.nvar = nvar;
newdata.p_id = p_id;
newdata.np = np;

% setup bounds
if isfield(optimdata,'lb')
    lb = optimdata.lb;
else
    lb = zeros(nvar,1);
end
if isfield(optimdata,'ub')
    ub = optimdata.ub;
else
    ub = zeros(nvar,1);
end
if isfield(optimdata,'eps_c')
    eps_c = optimdata.eps_c;
end
if isfield(optimdata,'eps_v')
    eps_v = optimdata.eps_v;
end
% set bounds
if type==1
    % set bounds - concentration
    lb(1:nc*npert) = 0; % xss_exp_v.*(1-eps_c);
    ub(1:nc*npert) = 10000; % xss_exp_v.*(1+eps_c);
    % set bounds - flux
    lb(nvar-nf*npert+1:nvar) = 0; % vss_exp_v*(1-eps_v);
    ub(nvar-nf*npert+1:nvar) = 10000; % vss_exp_v*(1+eps_v);
elseif type==2
    % set bounds - concentration
    lb(1:nc*npert) = xss_exp_v.*(1-eps_c);
    ub(1:nc*npert) = xss_exp_v.*(1+eps_c);
    % set bounds - flux
    lb(nc*npert+np+1:nc*npert+np+nf*npert) = vss_exp_v*(1-eps_v);
    ub(nc*npert+np+1:nc*npert+np+nf*npert) = vss_exp_v*(1+eps_v);
    % set bounds - noise
    lb(nvar-2+1:nvar) = 0;
    ub(nvar-2+1:nvar) = .5;
end
% set general bounds - parameter - specific bounds set below
lb(nc*npert+1:nc*npert+np) = .05*ones(np,1); 
ub(nc*npert+1:nc*npert+np) = 5*ones(np,1);

% setup fresh bounds specifically for each flux parameter dependent on
% flxid
if isfield(optimdata,'bounds')
    boundhs = optimdata.bounds;
end
if ~isempty(boundhs{flxid})
    boundfh = str2func(boundhs{flxid});
    [lb,ub] = boundfh(lb,ub,newdata);
end

% setup objectives and cons
if isfield(optimdata,'obj')
    objhs = optimdata.obj;
end
if isfield(optimdata,'nlcons')
    conshs = optimdata.nlcons;
end
if isfield(optimdata,'nlrhs')
    rhsval = optimdata.nlrhs;
else
    rhsval = [];
end
if isfield(optimdata,'nle')
    nles = optimdata.nle;
else
    nles = [];
end

newdata = rmfield(newdata,{'obj','nlcons','nlrhs','nle'});

% choose obj and cons
fhobj = str2func(objhs{flxid});
obj = @(x)fhobj(x,newdata.odep,newdata);
newdata.obj = obj;

fhcons = str2func(conshs{flxid});
nlcons = @(x)fhcons(x,newdata.odep,newdata);
newdata.nlcons = nlcons;

if ~isempty(rhsval{flxid})
    nlrhs = rhsval{flxid};
else
    if type==1
        nlrhs = zeros(npert,1);
    elseif type==2
        nlrhs = zeros(npert+2*nc*npert+2*npert,1);
    end
end
newdata.nlrhs = nlrhs;
if ~isempty(nles{flxid})
    nle = nles{flxid};
else
    if type==1
        nle = zeros(npert,1);
    elseif type==2
       nle = [zeros(npert,1);-ones(2*nc*npert+2*npert,1)];
    end
end
newdata.nle = nle;

% setup problem structure
prob =...
struct('obj',obj,'nlcons',nlcons,'nlrhs',nlrhs,'nle',nle,'lb',lb,'ub',ub);



    