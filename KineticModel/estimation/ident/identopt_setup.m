function prob_cas = identopt_setup(data,fixed_pvalue)
% p - all parameters j!=i optimized for
% ident_c - constant parameter identified
% p_useless - useless parameters not currently used in model
% acetate - fixed parametera acetate concentrations

if isfield(data,'casmodelfun')
    casmodelf = data.casmodelfun;
end
if isfield(data,'integratorfun')
    intfun = data.integratorfun;
end
if isfield(data,'xinit')
    xinit = data.xinit;
end
if isfield(data,'xexp')
    xexp = data.xexp;
end
if isfield(data,'tspan')
    tspan = data.tspan;
end
npts = length(tspan)-1;

% create CAS function with custom integrator of nlsq opt
intfun = str2func(intfun);
[xdyn_fun,~,~,p,ident_c,p_useless,acetate] = intfun(casmodelf,data);

% final symbolic expression to be used during optimization
x_sym =...
xdyn_fun(xinit,repmat(p,1,npts),fixed_pvalue,repmat(data.odep(14:16)',1,npts),.1);
% add initial value
x_sym = [casadi.DM(xinit) x_sym];
% choose only points present in experimental data
x_model_sym = x_sym(:,1:100:2001);

% create nlsqopt objective function
x_error = (xexp-x_model_sym);
obj = .5*dot(x_error,x_error);

prob_cas = struct('obj',obj,'x',p);
