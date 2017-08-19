function allfh = compare_vals(est_sol,exp_sol,data,opts,test_data_id)
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
% opt_xss(:,1) = [];
% opt_fss(:,1) = [];

% plot comparison
% use cell arrays to store data for each subplot
xss_plot = cell(data.nc,1);
xss_error = cell(data.nc,1);
fss_plot = cell(nf,1);
fss_error = cell(nf,1);

barc = struct();
for j = 1:data.nc
    xss_plot{j} = [exp_xss(j,:)' est_xss(j,:)'];
    xss_error{j} = [zeros(1,length(exp_xss(j,:)))' est_xerr(j,:)'];
end

for k = 1:nf
    fss_plot{k} = [exp_fss(k,:)' est_fss(k,:)'];
    fss_error{k} = [zeros(1,length(exp_fss(k,:)))' est_ferr(k,:)'];
end
hfc = figure;
ahc = zeros(data.nc,1);
for j = 1:data.nc
    ahc(j) = subplot(data.nc,1,j);
    set(ahc(j),'NextPlot','add');
    barc(j).h = bar(ahc(j),xss_plot{j});   
    [~,ylbl] = getKotteaxislabels(2,2,[1,j]);
    if strcmpi(version('-release'),'2014a')
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')
        
        ahc(j).YLabel.String = ylbl;    
        ahc(j).XTickLabel = rel_labels;  
    end
end

% collect data for error bars
% err_data = struct([]);
for j = 1:data.nc
    if strcmpi(version('-release'),'2014a')
        xpos = zeros(size(xss_plot{j},1),length(barc(j).h));
        ypos = zeros(size(xss_plot{j},1),length(barc(j).h));
        for ih = 1:length(barc(j).h)
            xpos_d = get(get(barc(j).h(ih),'children'),'xdata');
            ypos_d = get(get(barc(j).h(ih),'children'),'ydata');
            xpos(:,ih) = ((xpos_d(2,:)+xpos_d(3,:))/2)';
            ypos(:,ih) = ypos_d(2,:)';
        end
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')
        xpos = (cat(1,barc(j).h.XData)+cat(1,barc(j).h.XOffset))';
        ypos = xss_plot{j};
            
    end
    barc(j).xpos = xpos;
    barc(j).ypos = ypos;
end

% plot error bars
for j = 1:data.nc
    set(hfc,'CurrentAxes',ahc(j));
    errorbar(barc(j).xpos,barc(j).ypos,xss_error{j},'LineStyle','none');
end

% ahc.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');

% fluxes
hfv = figure;
ahf = zeros(nf,1);
barf = struct();
for k = 1:nf
    ahf(k) = subplot(nf/2,2,k);    
    set(ahf(k),'NextPlot','add');
    barf(k).h = bar(ahf(k),fss_plot{k});
    [~,ylbl] = getKotteaxislabels(2,1,[1,k]);
    if strcmpi(version('-release'),'2014a')
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')        
        ahf(k).YLabel.String = ylbl;  
        ahf(k).XTickLabel = rel_labels;
    end
end
% collect dsata for error bars
for k = 1:nf
    if strcmpi(version('-release'),'2014a')
        xpos = zeros(size(fss_plot{k},1),length(barf(k).h));
        ypos = zeros(size(fss_plot{k},1),length(barf(k).h));
        for ih = 1:length(barf(j).h)
            xpos_d = get(get(barf(k).h(ih),'children'),'xdata');
            ypos_d = get(get(barf(k).h(ih),'children'),'ydata');
            xpos(:,ih) = ((xpos_d(2,:)+xpos_d(3,:))/2)';
            ypos(:,ih) = ypos_d(2,:)';
        end
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')
        xpos = (cat(1,barf(k).h.XData)+cat(1,barf(j).h.XOffset))';
        ypos = fss_plot{j};
            
    end
    barf(k).xpos = xpos;
    barf(k).ypos = ypos;
end

% plot error bars
for k = 1:nf
    set(hfv,'CurrentAxes',ahf(k));
    errorbar(barf(k).xpos,barf(k).ypos,fss_error{k},'LineStyle','none');
end
% ahf.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Model Estimate');
allfh = [hfc;hfv];
