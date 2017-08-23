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
est_xss = cat(1,optsol.xss);
% est_fss = cat(2,est_data.fss);
est_par = cat(2,optsol.xpar);

% combine all experimental data
exp_xss = cat(2,exp_data.xss);
% exp_fss = cat(1,exp_data.fss);

wt_exp_xss = exp_xss(:,end); % wt exp data
wt_est_xss = est_xss(:,end); % wt est data

% re-do perturbations in exp_data to check for consistency of solutions
% initial value for these perturbations is taken from the wt model estimate
% create parameter and options structures for these perturbations
opts = struct([]);
newodep = repmat(data.odep,nval,1);
newodep(:,data.p_id) = est_par';
for ival = 1:nval    
    opts(ival).x0 = wt_est_xss(3*(ival-1)+1:3*ival);
    opts(ival).tspan = odeopts.tspan;
    opts(ival).solver_opts = odeopts.solver_opts;
    opts(ival).odep = newodep(ival,:);    
end
npert = size(exp_data,2);
[pt_val(1:npert).exp_pid] = exp_data(1:npert).exp_pid;
[pt_val(1:npert).exp_pval] = exp_data(1:npert).exp_pval;

parfor ival = 1:nval
    pss = ones(1,npert);
    % calculate fss for estimated xss
    odep = newodep(ival,:); 
    wt_est_fss = kotte_flux_noCAS(wt_est_xss(3*(ival-1)+1:3*ival),odep);
    
    new_odeopts = opts(ival);
    % do wt perturbation first to get ss (if it already is not @ ss)
    sol_wt = getperturbations(pt_val(end),@perturb_nonoise,new_odeopts);
    new_odeopts.x0 = sol_wt.xss;
    
    % recalculate perturbations for new parameters from new wt ss
    sol = getperturbations(pt_val(1:end-1),@perturb_nonoise,new_odeopts);
    close all

    % collect all solutions
%     cell_optsol{ival}.xss = [cat(2,sol.xss) wt_est_xss];
%     cell_optsol{ival}.fss = [cat(2,sol.fss) wt_est_fss];
    xss = [cat(2,sol.xss) sol_wt.xss];
    fss = [cat(2,sol.fss) sol_wt.fss];
    optsol(ival).xss_calc = xss;
    optsol(ival).fss_calc = fss;
    pss(xss(2,:)>xss(1,:)) = 2;
    optsol(ival).pss = pss;
end

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