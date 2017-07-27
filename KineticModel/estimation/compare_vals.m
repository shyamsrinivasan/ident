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
nf = size(wt_fss,1);
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
xss_plot = zeros(np+1,2*data.nc);
fss_plot = zeros(np+1,2*size(wt_fss,1));
for j = 0:data.nc-1
    xss_plot(:,2*j+1:2*(j+1)) =...
    [data.wt_xss(j+1) opt_xss(j+1);exp_xss(j+1,:)' est_xss(j+1,:)'];
end
for k = 0:nf-1
    fss_plot(:,2*k+1:2*(k+1)) =...
    [data.wt_fss(k+1) opt_fss(k+1);exp_fss(k+1,:)' est_fss(k+1,:)'];
end
figure
for j = 0:data.nc-1
    ahc = subplot(data.nc,1,j+1);
    bh = bar(ahc,xss_plot(:,2*j+1:2*(j+1)));
    [~,ylbl] = getKotteaxislabels(2,2,[1,j+1]);
    ahc.YLabel.String = ylbl;    
    ahc.XTickLabel = {'WT','P1','P2','P3'};
end
ahc.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');
% fluxes
figure
for k = 0:nf-1
    ahf = subplot(nf/2,2,k+1);
    bh = bar(ahf,fss_plot(:,2*k+1:2*(k+1)));
    [~,ylbl] = getKotteaxislabels(2,1,[1,k+1]);
    ahf.YLabel.String = ylbl;  
    ahf.XTickLabel = {'WT','P1','P2','P3'};
end
ahf.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');
