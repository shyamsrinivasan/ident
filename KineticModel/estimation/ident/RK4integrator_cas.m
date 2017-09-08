% integrate kotte model w/ custom RK4 integrator using mapaccum function of
% CasADi
function [xdyn_fun,ode,x,p_all,ident_c,p_useless,acetate] =...
        RK4integrator_cas(casfh,data)

if isfield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'npts')
    npts = data.npts;
end
if isfield(data,'ident_idx')
    ident_idx = data.ident_idx;
end

[ode,flux,~,~,x,p_all,ident_c,p_useless,acetate] = casfh(ident_idx,data.nc,data.nf,data.npert);
% [ode,~,~,fx,x,p] = kotte_CAS();

% RK4 integrator
dt = tspan(end)/npts;
k1 = ode(x,p_all,ident_c,p_useless,acetate);
k2 = ode(x+dt/2.0*k1,p_all,ident_c,p_useless,acetate);
k3 = ode(x+dt/2.0*k2,p_all,ident_c,p_useless,acetate);
k4 = ode(x+dt*k3,p_all,ident_c,p_useless,acetate);
xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
yfinal = flux(xfinal,p_all,ident_c,p_useless,acetate);

xstate_one_step =...
casadi.Function('xstate_one_step',{x,p_all,ident_c,p_useless,acetate},{xfinal,yfinal});

[xstate,ystate] = xstate_one_step(x,p_all,ident_c,p_useless,acetate);
xstate_onepoint =...
casadi.Function('xstate_onepoint',{x,p_all,ident_c,p_useless,acetate},{xstate,ystate});
xstate_onepoint = xstate_onepoint.expand();
xdyn_fun = xstate_onepoint.mapaccum('all_samples',npts);
