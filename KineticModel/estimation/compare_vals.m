function compare_vals(optsol,exp_sol,data,opts,exp_data_id)

exp_data_id = logical(exp_data_id);

% parse data to get optimal concentrations and parameters
[opt_xss,xpar] = parsesolvec(optsol,data);
opt_odep = data.odep;
opt_odep(data.p_id) = xpar;

% calculate new model fluxes
% before perturbation - wt fluxes
% assume wt (unperturbed) data is always included at the end
wt_exp_xss = exp_sol(end).xss;
wt_exp_fss = kotte_flux_noCAS(wt_exp_xss,data.odep);
% before perturbation - pertubed fluxes(using estimate concentrations)
opt_fss = kotte_flux_noCAS(opt_xss,opt_odep);

% unperturbed model fluxes when wt fluxes are not estimated
if ~exp_data_id(end)
    wt_est_xss = wt_exp_xss; % model wt xss = exp wt xss
    wt_est_fss = kotte_flux_noCAS(wt_est_xss,opt_odep);
else
    wt_est_xss = opt_xss(:,end);
    wt_est_fss = kotte_flux_noCAS(wt_est_xss,opt_odep);
end

% after perturbation - perturb parameters to ascertain fluxes
% np = length(find(exp_data_id));
pt_datasize = length(find(exp_data_id(1:end-1)));
nf = size(wt_exp_fss,1);
opts.odep = opt_odep;
opts.tspan = 0:.1:500;
% if np == size(exp_sol,2)
%     [pt_val(1:np-1).exp_pid] = exp_sol(1:end-1).exp_pid;
%     [pt_val(1:np-1).exp_pval] = exp_sol(1:end-1).exp_pval;
% else
    [pt_val(1:pt_datasize).exp_pid] = exp_sol(exp_data_id(1:end-1)).exp_pid;
    [pt_val(1:pt_datasize).exp_pval] = exp_sol(exp_data_id(1:end-1)).exp_pval;
% end
sol = getperturbations(pt_val,@perturb_nonoise,opts);
close all

% collect all data to plot [p#1 p#2....  wt]  
exp_xss = cat(2,exp_sol.xss);
exp_xss = exp_xss(:,exp_data_id);
% exp_xss = [exp_xss(:,end) exp_xss(:,1:end-1)];
exp_fss = cat(2,exp_sol.fss);
exp_fss = exp_fss(:,exp_data_id);
% exp_fss = [exp_fss(:,end) exp_fss(:,1:end-1)];

if ~exp_data_id(end) % wt not included in estimation
    est_xss = [cat(2,sol.xss) wt_est_xss];
    est_fss = [cat(2,sol.fss) wt_est_fss];    
else
    est_xss = [cat(2,sol.xss) opt_xss(:,end)];
    est_fss = [cat(2,sol.fss) opt_fss(:,end)];
end
nplot = pt_datasize+1;
% opt_xss(:,1) = [];
% opt_fss(:,1) = [];

% plot comparison
xss_plot = zeros(nplot,2*data.nc);
fss_plot = zeros(nplot,2*size(wt_exp_fss,1));
for j = 0:data.nc-1
    xss_plot(:,2*j+1:2*(j+1)) =...
    [exp_xss(j+1,:)' est_xss(j+1,:)'];
end
for k = 0:nf-1
    fss_plot(:,2*k+1:2*(k+1)) =...
    [exp_fss(k+1,:)' est_fss(k+1,:)'];
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
