% function to plot 2-d and 3-d diagrams comapring parameters
function hfp = compare_pars(est_data,plotype)
if nargin<2
    plotype=1;
end
if isempty(est_data)
    fprintf('No optimal solution found\n');
    hfp = [];
    return
end

hfp = figure;
xpar = est_data.xpar;
[npar,~] = size(xpar);
par = est_data.par;
perr = est_data.perr;

if plotype==1 % plot a bar graph of parameter value averages and error bars
    set(gca,'NextPlot','add');
    bh = bar(par);
    if strcmpi(version('-release'),'2014a')
        xpos_d = get(get(bh,'children'),'xdata');
        ypos_d = get(get(bh,'children'),'ydata');
        xerr_pos = ((xpos_d(2,:)+xpos_d(3,:))/2)';
        yerr_pos = ypos_d(2,:)';
    elseif strcmpi(version('-release'),'2014b') ||...
           strcmpi(version('-release'),'2017a') ||...
           strcmpi(version('-release'),'2017b')
        xerr_pos = cat(1,bh.XData)+cat(1,bh.XOffset);
        yerr_pos = par;
    end
    erh = errorbar(xerr_pos,yerr_pos,perr,'LineStyle','none');
elseif plotype==2
    if npar<=2 % 2-d plot
        line(xpar(1,:),xpar(2,:),'LineStyle','none',...
                                'Marker','.',...
                                'MarkerSize',10);
    elseif npar<=3 % 3-d plot
        line(xpar(1,:),xpar(2,:),xpar(3,:),'LineStyle','none',...
                                        'Marker','.',...
                                        'MarkerSize',10);
    end
end


