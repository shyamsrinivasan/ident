% integrate kotte model w/ custom RK4 integrator using mapaccum function of
% CasADi
function [xdyn_fun,ode,x,p] = RK4integrator_cas(casfh)

[ode,~,~,~,x,p] = casfh();
% [ode,~,~,fx,x,p] = kotte_CAS();

% RK4 integrator
dt = .1;
k1 = ode(x,p);
k2 = ode(x+dt/2.0*k1,p);
k3 = ode(x+dt/2.0*k2,p);
k4 = ode(x+dt*k3,p);

xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
xstate_one_step = casadi.Function('xstate_one_step',{x,p},{xfinal});

xstate = xstate_one_step(x,p);
xstate_onepoint = casadi.Function('xstate_onepoint',{x,p},{xstate});
xstate_onepoint = xstate_onepoint.expand();

xdyn_fun = xstate_onepoint.mapaccum('all_samples',3000);

% xdyn = all_samples(opts.x0,repmat(opts.odep(2:end)',1,3000))




