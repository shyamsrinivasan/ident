function proc_data =...
combine_results_v1(est_data,exp_data,val_data,data,odeopts)

if ~isfield(est_data,'xconc')
    est_data = parse_optsol(est_data,data);
end

% collect only optimal solutions
exitflags = cat(1,est_data.exitflag);
optsol = est_data(exitflags==1);
nval = size(optsol,2); % recalculate nval for optimal solutions only

% combine all estimated data
% est_xss = cat(1,optsol.xss);
% est_fss = cat(2,est_data.fss);
% est_par = cat(2,optsol.xpar);

% combine all experimental data
exp_xss = cat(2,exp_data.xss);
% exp_fss = cat(1,exp_data.fss);

wt_exp_xss = exp_xss(:,end); % wt exp data
% wt_est_xss = est_xss(:,end); % wt est data

% re-do perturbations in exp_data to check for consistency of solutions
% initial value for these perturbations is taken from the wt model estimate
% create parameter and options structures for these perturbations
optsol = reestimate_optsol(optsol,exp_data,data.odep,data.p_id,odeopts);

% collect all concentrations and fluxes
npar = length(data.p_id);

conc = cat(1,optsol.xss_calc);
flux = cat(1,optsol.fss_calc);
par = cat(2,optsol.xpar);
proc_data.calc_conc = conc;
proc_data.calc_flux = flux;
proc_data.par = par;
proc_data.pss = cat(1,optsol.pss);

% calculate averages and standard deviations based on mult
[avg,sigma2] = calc_avgstdev(nval,conc,data.nc,flux,data.nflx,par,npar);
proc_data.calc_xss = avg.avg_x;
proc_data.calc_xerr = sigma2.sigma2_x;
proc_data.calc_fss = avg.avg_f;
proc_data.calc_ferr = sigma2.sigma2_f;
proc_data.calc_par = avg.avg_p;
proc_data.calc_perr = sigma2.sigma2_p;

% separate data from different steady states
proc_data.calc_xss1
% proc_data.calc_xss2