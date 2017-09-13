function [prob_cas,data] = identopt_setup(data,fixed_pvalue)
% p - all parameters j!=i optimized for
% ident_c - constant parameter identified
% p_useless - useless parameters not currently used in model
% acetate - fixed parametera acetate concentrations

if isfield(data,'pname')
    pname = data.pname;
end
if isfield(data,'casmodelfun')
    casmodelf = data.casmodelfun;
end
if isfield(data,'integratorfun')
    intfun = data.integratorfun;
end
if isfield(data,'xinit')
    xinit = data.xinit;
end
if isfield(data,'yinit')
    yinit = data.yinit;
end    
if isfield(data,'xexp')
    xexp = data.xexp;    
end
if isfield(data,'yexp')
    yexp = data.yexp;
else
    yexp = xexp;
end
if isfield(data,'ynoise')
    % measurement noise/error => y ~ N(0,sigma^2);
    ynoise = data.ynoise;
else
    ynoise = ones(size(yexp,1),size(yexp,2));
%     y_noise(y_noise==1) = ;
end
if isfield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'freq')
    freq = data.freq;
end
npts = length(tspan)-1;
data.npts = npts;

if ~isfield(data,'ident_idx') && ~isempty(pname)    
    plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep',...
             'V4max','k1cat','V3max','V2max',...
             'vemax','KeFDP','ne','d','acetate'}; 
    ident_idx = find(strcmpi(plist,pname));    
    data.ident_idx = ident_idx;
elseif ~isfield(data,'ident_idx') && isempty(pname)    
    error('No parameter chosen for identifiability analysis');
end    

% create CAS function with custom integrator of nlsq opt
intfun = str2func(intfun);
[xdyn_fun,~,x,p_var,p_ident,p_fixed,acetate] = intfun(casmodelf,data);

% final symbolic expression to be used during optimization
[x_sym,y_sym] =...
xdyn_fun(xinit,repmat(p_var,1,npts),...
               fixed_pvalue,repmat(data.odep(6:9)',1,npts),...
               data.odep(17));

% x_sym = x_sym; % [sum(x_sym(1:3,:));sum(x_sym(4:6,:));sum(x_sym(7:9,:))];
% add initial value
x_sym = [casadi.DM(xinit) x_sym];
y_sym = [casadi.DM(yinit) y_sym];
% y_sym = [casadi.DM([sum(xinit(1:3,:));sum(xinit(4:6,:));sum(xinit(7:9,:))]) y_sym];
% choose only points present in experimental data
% x_model_sym = x_sym(:,freq);
% y_model_sym = [x_sym(1:2,freq);y_sym([1 3 4 5],freq)];
y_model_sym = y_sym([1 3 4 5],freq);
% create nlsqopt objective function
% x_error = (xexp-x_model_sym);
y_error = (yexp-y_model_sym)./ynoise;
obj = .5*dot(y_error,y_error);

objfun = casadi.Function('objfun',{p_var},{obj});

[lb,ub] = ident_bounds(length(p_var));

% jac = jacobian(obj,p);
jacfun = []; % casadi.Function('jacfun',{x,p,ident_c,p_useless,acetate},{jac});

% hess = hessian(obj,p);
hessfun = []; % casadi.Function('hessfun',{x,p,ident_c,p_useless,acetate},{hess});

prob_cas = struct('obj',obj,'x',x,'p',p_var,'ident_c',p_ident,....
                'p_useless',p_fixed,'acetate',acetate,'xdynfun',xdyn_fun,...
                'objfun',objfun,'grad',jacfun,'hess',hessfun,...
                'npts',npts,'xinit',xinit,'yinit',yinit,...
                'xexp',xexp,'yexp',yexp,...
                'lb',lb,'ub',ub);
