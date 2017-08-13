% integrate kotte model w/ custom RK4 integrator using mapaccum function of
% CasADi
function [xdyn_fun,ode,x,p_all,ident_c,p_useless,acetate] =...
        RK4integrator_cas(casfh,data)

if isfield(data,'pname')
    pname = data.pname;
end
if isfield(data,'tspan')
    tspan = data.tspan;
end
npts = length(tspan)-1;

plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
        'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA','acetate'}; 
if ~isempty(pname)
    ident_idx = find(strcmpi(plist,pname));
else
    ident_idx = [];
end    

[ode,~,~,~,x,p_all,ident_c,p_useless,acetate] = casfh(ident_idx);
% [ode,~,~,fx,x,p] = kotte_CAS();

% RK4 integrator
dt = tspan(end)/npts;
k1 = ode(x,p_all,ident_c,p_useless,acetate);
k2 = ode(x+dt/2.0*k1,p_all,ident_c,p_useless,acetate);
k3 = ode(x+dt/2.0*k2,p_all,ident_c,p_useless,acetate);
k4 = ode(x+dt*k3,p_all,ident_c,p_useless,acetate);

xfinal = x+dt/6.0*(k1+2*k2+2*k3+k4);
xstate_one_step =...
casadi.Function('xstate_one_step',{x,p_all,ident_c,p_useless,acetate},{xfinal});

xstate = xstate_one_step(x,p_all,ident_c,p_useless,acetate);
xstate_onepoint =...
casadi.Function('xstate_onepoint',{x,p_all,ident_c,p_useless,acetate},{xstate});
xstate_onepoint = xstate_onepoint.expand();

xdyn_fun = xstate_onepoint.mapaccum('all_samples',npts);

% xdyn = all_samples(opts.x0,repmat(opts.odep(2:end)',1,3000));




