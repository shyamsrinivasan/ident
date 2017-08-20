% scatter plot of experimental and estimated data
function compare_vals_scatter(est_sol,exp_sol,data,opts,test_data_id)
if isempty(est_sol)
    fprintf('No optimal solution found\n');
    allfh = [];
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
for j = 1:data.nc
    xss_plot{j} = [exp_xss(j,:)' est_xss(j,:)'];
    xss_error{j} = [zeros(1,length(exp_xss(j,:)))' est_xerr(j,:)'];
end

for k = 1:nf
    fss_plot{k} = [exp_fss(k,:)' est_fss(k,:)'];
    fss_error{k} = [zeros(1,length(exp_fss(k,:)))' est_ferr(k,:)'];
end

hfc = figure
