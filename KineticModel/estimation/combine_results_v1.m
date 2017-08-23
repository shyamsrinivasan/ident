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
% 
% exp_fss = cat(1,exp_data.fss);

% wt_exp_xss = exp_xss(:,end); % wt exp data
% wt_est_xss = est_xss(:,end); % wt est data

% re-do perturbations in exp_data to check for consistency of solutions
% initial value for these perturbations is taken from the wt model estimate
% create parameter and options structures for these perturbations
if nval>0
    optsol = reestimate_optsol(optsol,exp_data,data.odep,data.p_id,odeopts);
else
    proc_data = [];
    return
end

% collect all concentrations and fluxes
proc_data.opt_xss = cat(1,optsol.xss);
proc_data.opt_fss = cat(1,optsol.fss);
proc_data.opt_pss = cat(1,optsol.pss_opt);

npar = length(data.p_id);

conc = cat(1,optsol.xss_calc);
flux = cat(1,optsol.fss_calc);
par = cat(2,optsol.xpar);
proc_data.calc_xss = conc;
proc_data.calc_fss = flux;
proc_data.p = par;
proc_data.calc_pss = cat(1,optsol.pss_calc);

% calculate averages and standard deviations based on mult
[avg,sigma2] = calc_avgstdev(nval,conc,data.nc,flux,data.nflx,par,npar);
proc_data.calc_xss_avg = avg.avg_x;
proc_data.calc_xss_err = sigma2.sigma2_x;
proc_data.calc_fss_avg = avg.avg_f;
proc_data.calc_fss_err = sigma2.sigma2_f;
proc_data.calc_p_avg = avg.avg_p;
proc_data.calc_p_err = sigma2.sigma2_p;

% separate data from different steady states
% calculate ss of experimental data
exp_xss = cat(2,exp_data.xss);
exp_pss = ones(1,size(exp_data,2));
exp_pss(exp_xss(2,:)>exp_xss(1,:)) = 2;
exp_data.pss = exp_pss;

proc_data.calc_xss1
% proc_data.calc_xss2