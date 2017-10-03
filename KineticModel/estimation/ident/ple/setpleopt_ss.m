function [prob_cas,data] = setpleopt_ss(data,fixed_pvalue)
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
if isfield(data,'ynoise_var')
    % measurement noise/error => y ~ N(0,sigma^2);
    ynoise_var = data.ynoise_var;
end
if isfield(data,'input_data')
    input_data = data.input_data;
end
if isfield(data,'npts')
    npts = data.npts;
else
    npts = length(input_data);
end
if isfield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'freq')
    freq = data.freq;
end
% npts = length(tspan)-1;
% data.npts = npts;

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
[ode,flux,~,~,x,p_var,p_ident,p_fixed,acetate] =...
casmodelf(ident_idx,data.nc,data.nf,data.npert);

% RK4 integrator
dt = freq;
k1 = ode(x,p_var,p_ident,p_fixed,acetate);
k2 = ode(x+dt/2.0*k1,p_var,p_ident,p_fixed,acetate);
k3 = ode(x+dt/2.0*k2,p_var,p_ident,p_fixed,acetate);
k4 = ode(x+dt*k3,p_var,p_ident,p_fixed,acetate);
xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = flux(xfinal,p_var,p_ident,p_fixed,acetate);

% Create a function that simulates one step propagation in a sample
xstate_one_step =...
casadi.Function('xstate_one_step',...
                {x,p_var,p_ident,p_fixed,acetate},...
                {xfinal,yfinal});

% intfun = str2func(intfun);
% [xdyn_fun,~,x,p_var,p_ident,p_fixed,acetate] = intfun(casmodelf,data);
xstate = x;
for i=1:tspan(end)/freq
    [xstate,ystate] = xstate_one_step(xstate,p_var,p_ident,p_fixed,acetate);
end
% Create a function that simulates all step propagation on a sample
xstate_onepoint =...
casadi.Function('xstate_onepoint',...
                {x,p_var,p_ident,p_fixed,acetate},...
                {xstate,ystate});
xstate_onepoint = xstate_onepoint.expand();
xdyn_fun = xstate_onepoint.mapaccum('xdyn_fun',npts);

% final symbolic expression to be used during optimization
[x_sym,y_sym] =...
xdyn_fun(xinit,...
         repmat(p_var,1,npts),...
         fixed_pvalue,...
         repmat(data.odep(6:9)',1,npts),...
         input_data);

y_model_sym = y_sym([1 3 4 5],:);
% create nlsqopt objective function
% x_error = (xexp-x_model_sym);
y_error = ((yexp-y_model_sym).^2)./ynoise_var;
obj = .5*sum(sum(y_error));
objfun = casadi.Function('objfun',{p_var},{obj});

[lb,ub] = ident_bounds(length(p_var));

% jac = jacobian(obj,p);
jacfun = []; % casadi.Function('jacfun',{x,p,ident_c,p_useless,acetate},{jac});

% hess = hessian(obj,p);
hessfun = []; % casadi.Function('hessfun',{x,p,ident_c,p_useless,acetate},{hess});

prob_cas = struct('obj',obj,'p',p_var,'ident_c',p_ident,....
                'p_useless',p_fixed,'acetate',acetate,'xdynfun',xdyn_fun,...
                'objfun',objfun,'grad',jacfun,'hess',hessfun,...
                'npts',npts,'xinit',xinit,'yinit',yinit,...
                'xexp',xexp,'yexp',yexp,...
                'lb',lb,'ub',ub);
