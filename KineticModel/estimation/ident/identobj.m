% using IPOPT outise casadi
function objfun = identobj(xvar,ident_c,data)
if isfield(data,'xinit')
    xinit = data.xinit;
end
if isfield(data,'ident_idx')
    idx = data.ident_idx;
else
    idx = 13;
end
% if isfield(data,'odep')
%     odep = data.odep;
% end
if isfield(data,'tspan')
    tspan = data.tspan;
end
if isfield(data,'nc')
    nc = data.nc;
end
if isfield(data,'npert')
    npert = data.npert;
end
if isfield(data,'nlaefh')
    nlaefh = data.nlaefh;
else
    nlaefh = @kotteCAS_pert;
end
if isfield(data,'odefh')
    casfh = data.odefh;
else
    casfh = @kotteCAS_pert;
end
% p_useless = data.odep(14:16)';
p_pert = reshape(data.p_pert,3*npert,1);
p_pert(14-idx) = ident_c;

ac = data.odep(end);
p_all = [xvar(1:idx-1);ident_c;xvar(idx:end);ac;p_pert];
% dynamic case
options = struct('grid',tspan,'output_t0',1,'print_stats',1);
ymodel = solve_model_ode(casfh,nc,npert,xinit,p_all,options);
% ss case
ymodel = solve_model_nlae(nlaefh,nc,npert,xinit,p_all);
yerror = yexpt-ymodel;
objfun = .5*dot(yerror,yerror);