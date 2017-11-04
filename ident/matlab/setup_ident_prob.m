% setup PLE problem for identifiability
function [prob,newdata] = setup_ident_prob(data)

if isfield(data,'pname')
    % index of parameter for identifiability analysis
    pname = data.pname;
else
    error('No parameters to estimate');
end

% all parameters available (ordered)
plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA','acetate'}; 
if ~isempty(pname)
    idx = find(strcmpi(plist,pname));
else
    idx = [];
end

% determine nvar for optimization
% nvar = # parameters to be optimized
if ~isempty(idx)
    nvar = length(plist)-1;    
end
newdata = data;
newdata.idx = idx;   
newdata.np = length(plist);
newdata.nvar = nvar;
newdata.varid = setdiff(1:length(plist),idx);

% setup parameter estimation bounds - specific bounds set below
if isfield(data,'lb')
    lb = data.lb;
else
    lb = zeros(nvar,1);
end
if isfield(data,'ub')
    ub = data.ub;
else
    ub = ones(nvar,1);
end

% setup fresh bounds specifically for each flux parameter 
if isfield(data,'boundfh')
    boundfh = data.boundfh;
end
if ~isempty(boundfh)
    boundfh = str2func(boundfh);
    [lb,ub] = boundfh(lb,ub,newdata);
end

% model data fun handle
if isfield(data,'modelf')
    modelfh = data.modelf;
else
    modelfh = [];
end
if ~isempty(modelfh)
    modelfh = str2func(modelfh);
    modelf = @(p)modelfh(p,newdata);
end
% model sensitivity fun handle
if isfield(data,'modelsensf')
    modelsensfh = data.modelsensf;
else
    modelsensfh = [];
end
if ~isempty(modelsensfh)
    modelsensfh = str2func(modelsensfh);
    modelsensf = @(p)modelsensfh(p,newdata);
end

newdata.modelf = modelf;
newdata.modelsensf = modelsensf;
newdata = rmfield(newdata,{'pname','boundfh','objf','gradobjfh'});

% setup objectives
if isfield(data,'objf')
    objfh = data.objf;
else
    objfh = [];
end
if ~isempty(objfh)
    objfh = str2func(objfh);
    obj = @(p)objfh(p,newdata);
end

% setup objective gradients
if isfield(data,'gradobjfh')
    gradobjfh = data.gradobjfh;
else
    gradobjfh = [];
end
if ~isempty(gradobjfh)
    gradobjfh = str2func(gradobjfh);
    gradobj = @(p)gradobjfh(p,newdata);
end

% setup cons
if isfield(data,'nlconsf')
    nlcons = data.nlconsf;
    newdata = rmfield(newdata,'nlconsf');
else
    nlcons = [];
end
if isfield(data,'nlrhs')
    rhsval = data.nlrhs;
    newdata = rmfield(newdata,'nlrhs');
else
    rhsval = [];
end
if isfield(data,'nle')
    nle = data.nle;
    newdata = rmfield(newdata,'nle');
else
    nle = [];
end

prob =...
struct('obj',obj,'gradobj',gradobj,...
        'nlcons',nlcons,'nlrhs',rhsval,'nle',nle,...
        'lb',lb,'ub',ub);
