% scatter plot of experimental and estimated data
% plot all data in proc_data struct
function compare_vals_scatter(exp_sol,nexp,opt_xss,calc_xss,data,fc)
if nargin<8
    fc = 2; % concentrations plots
end

npert = data.npert;
% experimental data
if fc==2
    nvar = data.nc;
    nval = size(opt_xss,1)/nvar;
    exp_xss = extract_expdata(exp_sol,1,[],[],nvar,npert);
    % optimal data curation into different conc/var
    xss = cell(nvar,1);
    for j = 1:nvar    
        xss{j,1} = exp_xss(j,:);
        xss{j,2} = opt_xss(j:nvar:nvar*nval,:);
        xss{j,3} = calc_xss(j:nvar:nvar*nval,:);         
    end
elseif fc==1
    nvar = data.nflx;
    nval = size(opt_xss,1)/nvar;
    [~,exp_fss] = extract_expdata(exp_sol,1,[],nvar,[],npert);
    xss = cell(nvar,1);
    for j = 1:nvar    
        xss{j,1} = exp_fss(j:nvar:nvar*nval,:);
        xss{j,2} = opt_xss(j:nvar:nvar*nval,:);
        xss{j,3} = calc_xss(j:nvar:nvar*nval,:);        
    end
end

% plot scatter of xss{j,1} vs xss{j,2} and xss{j,3}
plotdots(nvar,xss(:,1),xss(:,2),xss(:,3),fc);
    

function [hf,dotfs] = plotdots(nplot,ss_x,ss_y1,ss_y2,fc)
if nargin<5
    fc = 2;
end

hf = figure;
dotfs = struct();
for j = 1:nplot
    if fc==2
        dotfs(j).ah = subplot(nplot,1,j);        
    elseif fc==1
        dotfs(j).ah = subplot(nplot/2,2,j);            
    end
    set(dotfs(j).ah,'NextPlot','add');
%     line(dotfs(j).ah,ss_x{j},ss_x{j},'Color','k','LineWidth',2);
    if ~isempty(ss_y1{j})        
        dotfs(j).h = line(dotfs(j).ah,ss_x{j},ss_y1{j},...
                          'LineStyle','none','Marker','.',...
                          'MarkerSize',18);
    end
    if ~isempty(ss_y2{j})
        dotfs(j).h = line(dotfs(j).ah,ss_x{j},ss_y2{j},...
                          'LineStyle','none','Marker','.',...
                          'MarkerSize',18,'Color','r');
    end
    [~,ylbl] = getKotteaxislabels(2,fc,[1,j]);
    if strcmpi(version('-release'),'2014a')
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')

        dotfs(j).ah.YLabel.String = ylbl{1};    
%         ah(j).XTickLabel = rel_labels;  
    end    
end




