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
close all

% collect all data to plot 
exp_xss = cat(2,exp_sol.xss);
exp_fss = cat(2,exp_sol.fss);

est_xss = cat(2,sol.xss);
est_fss = cat(2,sol.fss);

% plot comparison
x1ss = [exp_xss(1,:);est_xss(1,:)]';
x2ss = [exp_xss(2,:);est_xss(2,:)]';
x3ss = [exp_xss(3,:);est_xss(3,:)]';
figure
subplot(311)
bar(x1ss);
legend('Noisy Data','Model Estimate');
[~,ylbl] = getKotteaxislabels(2,2,[1,1]);
ylabel(ylbl);
subplot(312)
bar(x2ss);
[~,ylbl] = getKotteaxislabels(2,2,[1,2]);
ylabel(ylbl);
subplot(313)
bar(x3ss);
[~,ylbl] = getKotteaxislabels(2,2,[1,3]);
ylabel(ylbl);

f1ss = [exp_fss(1,:);est_fss(1,:)]';
f2ss = [exp_fss(2,:);est_fss(2,:)]';
f3ss = [exp_fss(3,:);est_fss(3,:)]';
f4ss = [exp_fss(4,:);est_fss(4,:)]';
f5ss = [exp_fss(5,:);est_fss(5,:)]';
f6ss = [exp_fss(6,:);est_fss(6,:)]';
% fluxes
figure
subplot(331)
bh = bar(f1ss);
legend('Noisy Data','Model Estimate');
[~,ylbl] = getKotteaxislabels(2,1,[1,1]);
ylabel(ylbl);
subplot(322)
bar(f2ss);
[~,ylbl] = getKotteaxislabels(2,1,[1,2]);
ylabel(ylbl);
subplot(323)
bar(f3ss);
[~,ylbl] = getKotteaxislabels(2,1,[1,3]);
ylabel(ylbl);
subplot(324)
bar(f4ss);
[~,ylbl] = getKotteaxislabels(2,1,[1,4]);
ylabel(ylbl);
subplot(325)
bar(f5ss);
[~,ylbl] = getKotteaxislabels(2,1,[1,5]);
ylabel(ylbl);
subplot(326)
bar(f6ss);
[~,ylbl] = getKotteaxislabels(2,1,[1,6]);
ylabel(ylbl);
% figure
% bar(est_fss);

