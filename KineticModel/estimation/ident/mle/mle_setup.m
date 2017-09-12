function prob = mle_setup(data)

if isfield(data,'casmodelfun')
    casfh = data.casmodelfun;
end
if isfield(data,'xinit')
    xinit = data.xinit;
end
if isfield(data,'yinit')
    yinit = data.yinit;
end 
if isfield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'freq')
    freq = data.freq;
end
if isfield(data,'xexp')
    xexp = data.xexp;    
end
if isfield(data,'yexp')
    yexp = data.yexp;
end
if isfield(data,'ynoise')
    % measurement noise/error => y ~ N(0,sigma^2);
    ynoise = data.ynoise;
else
    ynoise = ones(size(yexp,1),size(yexp,2));
%     y_noise(y_noise==1) = .001;
end
npts = length(tspan)-1;

% create CAS function with custom integrator of nlsq opt
[ode,flux,~,~,x,p_var,p_other,acetate] = casfh(data.nc);

% RK4 integrator
dt = tspan(end)/npts;
k1 = ode(x,p_var,p_other,acetate);
k2 = ode(x+dt/2.0*k1,p_var,p_other,acetate);
k3 = ode(x+dt/2.0*k2,p_var,p_other,acetate);
k4 = ode(x+dt*k3,p_var,p_other,acetate);
xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = flux(xfinal,p_var,p_other,acetate);

xstate_one_step =...
casadi.Function('xstate_one_step',{x,p_var,p_other,acetate},{xfinal,yfinal});

[xstate,ystate] = xstate_one_step(x,p_var,p_other,acetate);
xstate_onepoint =...
casadi.Function('xstate_onepoint',{x,p_var,p_other,acetate},{xstate,ystate});
xstate_onepoint = xstate_onepoint.expand();
xdyn_fun = xstate_onepoint.mapaccum('all_samples',npts);

% intfun = str2func(intfun);
% [xdyn_fun,~,x,p,ident_c,p_useless,acetate] = intfun(casmodelf,data);
% p_all = [p(1:ident_idx-1);ident_c;p(ident_idx:end)]; % ;
p_fixed = [p_other;acetate];

% final symbolic expression to be used during optimization
% [x_sym,y_sym] =...
% xdyn_fun(xinit,repmat(p,1,npts),repmat(data.odep(14:16)',1,npts),.1);
[x_sym,y_sym] =...
xdyn_fun(xinit,repmat(p_var,1,npts),repmat(data.odep(6:9)',1,npts),data.odep(17));

% add initial value
x_sym = [casadi.DM(xinit) x_sym];
y_sym = [casadi.DM(yinit) y_sym];

y_model_sym = y_sym([1 3 4 5],freq);
y_error = (yexp-y_model_sym)./ynoise;
obj = .5*dot(y_error,y_error);

% define obj as function to be used outside casadi implementation of ipopt
objfun = casadi.Function('objfun',{p_var},{obj});

% fixed parameters in p_fixed = {'vemax','KeFDP','ne','d','acetate'}
fixed_p = [data.odep(6:9)';data.odep(17)]; 
% fixed_p = data.odep(17); 
objfh = @(p)full(objfun(p));

% input for jac fun is same as objfun (function whose jacobian is calculated)
gradfun = jacobian(objfun); % this generates a function for jacobian
gradfh = @(p)full(gradfun(p));

% bounds for mle estimate of all 13 parameters for given input (acetate)
[lb,ub] = ident_bounds_mle(length(p_var));
hessfun = [];

prob = struct('objfun',objfh,'grad',gradfh,'hess',hessfun,...
            'npts',npts,'xinit',xinit,'yinit',yinit,...
            'xexp',xexp,'yexp',yexp,...
            'lb',lb,'ub',ub);