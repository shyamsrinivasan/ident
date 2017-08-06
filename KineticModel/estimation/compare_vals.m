function allfh = compare_vals(est_sol,exp_sol,data,opts,test_data_id)
if isempty(est_sol)
    fprintf('No optimal solution found\n');
    return
end

test_data_id = logical(test_data_id);
all_labels = {'P1','P2','P3','WT'};

% calculate new model fluxes
% before perturbation - wt fluxes
% assume wt (unperturbed) data is always included at the end
wt_exp_xss = exp_sol(end).xss;
wt_exp_fss = kotte_flux_noCAS(wt_exp_xss,data.odep);

% after perturbation - perturb parameters to ascertain fluxes
pt_datasize = length(find(test_data_id(1:end-1)));
nf = size(wt_exp_fss,1);

exp_xss = cat(2,exp_sol.xss);
exp_fss = cat(2,exp_sol.fss);

% collect all data to plot [p#1 p#2....  wt]  
% estimated
est_xss = est_sol.xss;
est_xerr = est_sol.xerr;
est_fss = est_sol.fss;
est_ferr = est_sol.ferr;

if ~test_data_id(end) % wt not included in estimation
    rel_labels = all_labels(test_data_id);
else
    % experimental
    exp_xss = exp_xss(:,test_data_id);
    exp_fss = exp_fss(:,test_data_id);
    % axis labels
    rel_labels = [all_labels(test_data_id) all_labels(end)];
end
nplot = pt_datasize+1;
% opt_xss(:,1) = [];
% opt_fss(:,1) = [];

% plot comparison
xss_plot = zeros(nplot,2*data.nc);
xss_error = zeros(nplot,2*data.nc);
fss_plot = zeros(nplot,2*nf);
fss_error = zeros(nplot,2*nf);
for j = 0:data.nc-1
    xss_plot(:,2*j+1:2*(j+1)) =...
    [exp_xss(j+1,:)' est_xss(j+1,:)'];
    xss_error(:,2*j+1:2*(j+1)) =...
    [zeros(1,length(exp_xss(j+1,:)))' est_xerr(j+1,:)'];
end
for k = 0:nf-1
    fss_plot(:,2*k+1:2*(k+1)) =...
    [exp_fss(k+1,:)' est_fss(k+1,:)'];
    fss_error(:,2*k+1:2*(k+1)) =...
    [zeros(1,length(exp_fss(k+1,:)))' est_ferr(k+1,:)'];
end
hfc = figure;
for j = 0:data.nc-1
    ahc = subplot(data.nc,1,j+1);
    set(ahc,'NextPlot','add');
    bh = bar(ahc,xss_plot(:,2*j+1:2*(j+1)));    
    xerr_pos = cat(1,bh.XData)+cat(1,bh.XOffset);
    yerr_pos = xss_plot(:,2*j+1:2*(j+1));    
    erh = errorbar(xerr_pos',yerr_pos,...
                             xss_error(:,2*j+1:2*(j+1)),'LineStyle','none');
               
    [~,ylbl] = getKotteaxislabels(2,2,[1,j+1]);
    ahc.YLabel.String = ylbl;    
%     ahc.XTickLabel = rel_labels;    
end
ahc.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');
% fluxes
hfv = figure;
for k = 0:nf-1
    ahf = subplot(nf/2,2,k+1);
    set(ahf,'NextPlot','add');
    bh = bar(ahf,fss_plot(:,2*k+1:2*(k+1)));
    xerr_pos = cat(1,bh.XData)+cat(1,bh.XOffset);
    yerr_pos = fss_plot(:,2*k+1:2*(k+1));    
    erh = errorbar(xerr_pos',yerr_pos,...
                             fss_error(:,2*k+1:2*(k+1)),'LineStyle','none');
    [~,ylbl] = getKotteaxislabels(2,1,[1,k+1]);
    ahf.YLabel.String = ylbl;  
%     ahf.XTickLabel = rel_labels;
    clear bh
end
ahf.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');
allfh = [hfc;hfv];
