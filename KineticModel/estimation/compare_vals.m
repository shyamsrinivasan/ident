function [allfh,barh] = compare_vals(est_sol,exp_sol,val_data,data,type)
if nargin<5
    type = 1;
end
if isempty(est_sol)
    fprintf('No optimal solution found\n');
    allfh = [];
    barh = [];
    return
end

% get est data from fields in est_sol
% optimal sol
% if isfield(est_sol,'opt_xss')
%     opt_xss = est_sol.opt_xss;
% end
% if isfield(est_sol,'opt_fss')
%     opt_fss = est_sol.opt_fss;
% end
if isfield(est_sol,'opt_xss_avg')
    opt_xss_avg = est_sol.opt_xss_avg;
else
    opt_xss_avg = [];
end
if isfield(est_sol,'opt_xss_err')
    opt_xss_err = est_sol.opt_xss_err;
else
    opt_xss_err = [];
end
if isfield(est_sol,'opt_fss_avg')
    opt_fss_avg = est_sol.opt_fss_avg;
else
    opt_fss_avg = [];
end
if isfield(est_sol,'opt_fss_err')
    opt_fss_err = est_sol.opt_fss_err;
else
    opt_fss_err = [];
end

% recalculated values
% if isfield(est_sol,'calc_xss')
%     calc_xss = est_sol.calc_xss;
% end
% if isfield(est_sol,'calc_fss')
%     calc_fss = est_sol.calc_fss;
% end
if isfield(est_sol,'calc_xss_avg')
    calc_xss_avg = est_sol.calc_xss_avg;
end
if isfield(est_sol,'calc_xss_err')
    calc_xss_err = est_sol.calc_xss_err;
end
if isfield(est_sol,'calc_fss_avg')
    calc_fss_avg = est_sol.calc_fss_avg;
end
if isfield(est_sol,'calc_fss_err')
    calc_fss_err = est_sol.calc_fss_err;
end

% experimental data
if isfield(exp_sol,'xss')
    xss = cat(2,exp_sol.xss);
end
if isfield(exp_sol,'fss')
    fss = cat(2,exp_sol.fss);
end

if type==1 % bar chart
    % plot concentrations
    xss_bar_data =...
    barplotdata(data.nc,xss,[],opt_xss_avg,opt_xss_err,calc_xss_avg,calc_xss_err);
    xss_plot = xss_bar_data.ss;
    xss_error = xss_bar_data.err;
    [hfc,barc] = plotbars(data.nc,xss_plot,xss_error);
    
    % plot fluxes
    fss_bar_data =...
    barplotdata(data.nflx,fss,[],...
                [opt_fss_avg;zeros(data.nflx-1,length(opt_fss_avg))],...
                [opt_fss_err;zeros(data.nflx-1,length(opt_fss_err))],...
                calc_fss_avg,calc_fss_err);
    fss_plot = fss_bar_data.ss;
    fss_error = fss_bar_data.err;
    [hfv,barv] = plotbars(data.nflx,fss_plot,fss_error,1);
    
elseif type==2 % scatter plot
end
    
% ahc.XLabel.String = 'WT and Perturbations';
legend('Noisy Data','Optimal Estimate','Model SS');
allfh = [hfc;hfv];
barh = {barc;barv};

function [hf,barfs] = plotbars(nplot,ss_plot,error_plot,fc)
if nargin<4
    fc = 2;
end

hf = figure;
barfs = struct();
for j = 1:nplot
    if fc==2
        barfs(j).ah = subplot(nplot,1,j);        
    elseif fc==1
        barfs(j).ah = subplot(nplot/2,2,j);            
    end
    set(barfs(j).ah,'NextPlot','add');
    barfs(j).h = bar(barfs(j).ah,ss_plot{j});   
    [~,ylbl] = getKotteaxislabels(2,fc,[1,j]);
    if strcmpi(version('-release'),'2014a')
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')

        barfs(j).ah.YLabel.String = ylbl{1};    
%         ah(j).XTickLabel = rel_labels;  
    end
end

% collect data for error bars
% err_data = struct([]);
for j = 1:nplot
    if strcmpi(version('-release'),'2014a')
        xpos = zeros(size(ss_plot{j},1),length(barfs(j).h));
        ypos = zeros(size(ss_plot{j},1),length(barfs(j).h));
        for ih = 1:length(barfs(j).h)
            xpos_d = get(get(barfs(j).h(ih),'children'),'xdata');
            ypos_d = get(get(barfs(j).h(ih),'children'),'ydata');
            xpos(:,ih) = ((xpos_d(2,:)+xpos_d(3,:))/2)';
            ypos(:,ih) = ypos_d(2,:)';
        end
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')
        xpos = (cat(1,barfs(j).h.XData)+cat(1,barfs(j).h.XOffset))';
        ypos = ss_plot{j};

    end
    barfs(j).xpos = xpos;
    barfs(j).ypos = ypos;
end

% plot error bars
for j = 1:nplot
    set(hf,'CurrentAxes',barfs(j).ah);
    errorbar(barfs(j).xpos,barfs(j).ypos,error_plot{j},'LineStyle','none');
end
