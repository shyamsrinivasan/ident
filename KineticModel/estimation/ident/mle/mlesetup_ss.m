% steady state data for MLE
% setup function for ss MLE 
% uses casadi CAS system to build models 
function prob = mlesetup_ss(data)

if isfield(data,'casmodelfun')
    casfh = data.casmodelfun;
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
end
if isfield(data,'xnoise_var')
    xnoise_var = data.xnoise_var;
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

% create CAS function with custom integrator of nlsq opt
[ode,flux,~,~,x,p_var,p_other,acetate] = casfh(data.nc);

% RK4 integrator
dt = freq;
k1 = ode(x,p_var,p_other,acetate);
k2 = ode(x+dt/2.0*k1,p_var,p_other,acetate);
k3 = ode(x+dt/2.0*k2,p_var,p_other,acetate);
k4 = ode(x+dt*k3,p_var,p_other,acetate);
xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = flux(xfinal,p_var,p_other,acetate);

% Create a function that simulates one step propagation in a sample
xstate_one_step =...
casadi.Function('xstate_one_step',{x,p_var,p_other,acetate},{xfinal,yfinal});
xstate = x;
for i=1:tspan(end)/freq
    [xstate,ystate] = xstate_one_step(xstate,p_var,p_other,acetate);
end

% Create a function that simulates all step propagation on a sample
xstate_onepoint =...
casadi.Function('xstate_onepoint',{x,p_var,p_other,acetate},{xstate,ystate});
xstate_onepoint = xstate_onepoint.expand();
xdyn_fun = xstate_onepoint.mapaccum('xdyn_fun',npts);
% final symbolic expression to be used during optimization
[x_sym,y_sym] =...
xdyn_fun(xinit,repmat(p_var,1,npts),repmat(data.odep(6:9)',1,npts),input_data);
if isfield(data,'wcon') && data.wcon
    % add mets concentration 
    y_model_sym = [y_sym([1 3 4 5],:);x_sym([1 2],:)]; 
    y_error = (([yexp;xexp]-y_model_sym).^2)./ynoise_var;
else
    y_model_sym = y_sym([1 3 4 5],:);
    y_error = ((yexp-y_model_sym).^2)./ynoise_var;
end

obj = .5*sum(sum(y_error));
objfun = casadi.Function('objfun',{p_var},{obj});

% bounds for mle estimate of all 13 parameters for given input (acetate)
[lb,ub] = ident_bounds_mle(length(p_var));
gradfh = [];
hessfun = [];

prob = struct('obj',obj,'objfun',objfun,'p',p_var,'grad',gradfh,'hess',hessfun,...
            'npts',npts,'xinit',xinit,'yinit',yinit,...
            'xexp',xexp,'yexp',yexp,...
            'lb',lb,'ub',ub);
