function compare_vals(optsol,exp_sol,data,opts)

% parse data to get optimal concentrations and parameters
[opt_xss,xpar] = parsesolvec(optsol,data);
opt_odep = data.odep;
opt_odep(data.p_id) = xpar;

% calculate new model fluxes
% before perturbation
wt_fss = kotte_flux_noCAS(data.wt_xss,data.odep);
opt_fss = kotte_flux_noCAS(opt_xss,opt_odep);

% after perturbation - perturb parameters to ascertain fluxes
np = size(exp_sol,2);
opts.odep = opt_odep;
opts.tspan = 0:.1:500;
[pt_val(1:np).exp_pid] = exp_sol.exp_pid;
[pt_val(1:np).exp_pval] = exp_sol.exp_pval;
sol = getperturbations(pt_val,@perturb_nonoise,opts);

% collect all data to plot 
exp_xss = cat(2,exp_sol.xss);
exp_fss = cat(2,exp_sol.fss);

est_xss = cat(2,sol.xss);
est_fss = cat(2,sol.xss);

% plot comparison


