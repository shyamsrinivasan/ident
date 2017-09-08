% identopt_setup for using IPOPT outside of casadi
function ymodel = solve_model_ode(casfh,nc,npert,xinit,odep,options)

[ode,~,oderhs,~,xstate,p_all] = casfh(nc,npert);

npts = length(options.grid)-1;
% RK4 integrator
dt = options.grid(end)/npts;
k1 = ode(xstate,p_all);
k2 = ode(xstate+dt/2.0*k1,p_all);
k3 = ode(xstate+dt/2.0*k2,p_all);
k4 = ode(xstate+dt*k3,p_all);
xfinal = xstate+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = xfinal;

xstate_one_step =...
casadi.Function('xstate_one_step',{xstate,p_all},{yfinal});
xstate_one_step_call = xstate_one_step(xstate,p_all);
xstate_onepoint =...
casadi.Function('xstate_onepoint',{xstate,p_all},{xstate_one_step_call});
xstate_onepoint = xstate_onepoint.expand();
xdyn_fun = xstate_onepoint.mapaccum('all_samples',npts);

% xstate_sym = xdyn_fun(xinit,repmat(p_all,1,npts));
ymodel = [xinit full(xdyn_fun(xinit,repmat(odep,1,npts)))];

% ode_struct = struct('x',xstate,'ode',oderhs,'p',p_all);
% int = casadi.integrator('int','cvodes',ode_struct,options);
% sol = int('x0',xinit,'p',odep);
% ymodel = full(sol.xf);


