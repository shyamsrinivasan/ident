% identopt_setup for using IPOPT outside of casadi
function ymodel = solve_model_ode(casfh,nc,npert,xinit,odep,options)

[~,~,oderhs,~,xstate,p_all] = casfh(nc,npert);

ode_struct = struct('x',xstate,'ode',oderhs,'p',p_all);
int = casadi.integrator('int','cvodes',ode_struct,options);
sol = int('x0',xinit,'p',odep);
ymodel = full(sol.xf);


