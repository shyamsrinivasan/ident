function plotflux_scatter(model,flux,fluxind,hsubfig)
if nargin < 3
    nflux = length(model.Vind)+length(model.Vexind);
else
    nflux = length(fluxind);
end
if nargin < 4
    hsubfig = zeros(nflux,1);
end
nsmpl = size(flux,2);
fluxt = flux';
if rem(nflux,2)==0
    n = nflux/2;
else
    n = (nflux+1)/2;
end
%Randomly choose colors
load('C:\Users\shyam\Documents\Courses\CHE 1125 Project\Colors.mat');
Clrnd = floor(1 + (76-1)*rand(nflux,1));
ColorSpec = cell(length(Clrnd),1);
for j = 1:length(Clrnd)
    try 
        ColorSpec{j} = rgb(Colors{Clrnd(j)});
    catch
        ColorSpec{j} = rgb('Black');
    end
end

% hline = zeros(nflux,1);
%Scatter Plot all Fluxes
for ivar = 1:nflux
    if isfield(model,'Regulators');
        tfp = find(strcmpi(fluxind{ivar},model.Regulators));
        y_label1 = sprintf('Flux %s',model.Regulators{tfp});
    else
        tfp = find(strcmpi(fluxind{ivar},model.Enzyme));
        y_label1 = sprintf('Flux %s',model.Enzyme{tfp});
    end
    if any(tfp)
        if isempty(findobj('type','figure','Name','Steady State Fluxes'))
            hfig = figure('Name','Steady State Fluxes'); 
            figure(hfig);
        else
            hfig = findobj('type','figure','Name','Steady State Fluxes');
            figure(hfig);
        end
%         if tfp ~= model.Vupind
            if hsubfig(ivar) ~= 0
                hca = findobj(hsubfig(ivar),'type','axes');
                set(hfig,'CurrentAxes',hca);  
            else%subplot is unassigned
                hsubfig(ivar) = subplot(2,n,ivar);       
                hca = gca;
                %Make sure more plots can be added at end of loop
                set(hca,'NextPlot','add');     
                set(hsubfig(ivar),'NextPlot','add');
            end
%             hsubfig(tfp-1) = subplot(nflux,1,ivar);
            fluxd = fluxt(:,tfp);%flux of ivar
%             fluxd_p = fluxt(:,prxnid);
            X = zeros(length(fluxd),1);
            initf = fluxd(1);
            fluxd(1) = [];
            X(1) = 1;
            X(2:length(fluxd)+1) = 2;
            hline = line(X,[initf;fluxd]);
            set(hline,'LineStyle','none',...
                         'Marker','o',...
                         'MarkerSize',5,...
                         'MarkerFaceColor',ColorSpec{ivar},...
                         'MarkerEdgeColor',ColorSpec{ivar});

            set(get(gca,'YLabel'),'String',y_label1);  
            set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
            set(get(gca,'YLabel'),'FontSize',12);  

            xlabel = sprintf('Steady States');
            set(get(gca,'XLabel'),'String',xlabel);  
            set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
            set(get(gca,'XLabel'),'FontSize',12);
            Ydata = get(hline,'YData');
            
            maxY = max(Ydata);
            minY = min(Ydata);
            [minY,maxY] = fixaxisbounds(minY,maxY);

            set(gca,'YLim',[minY maxY]);
            AP.XColor = [.1 .1 .1];
            AP.YColor = [.1 .1 .1];
            AP.XLim = [0 max(X)+1];
            AP.XTick = 0:1:max(X);
            set(hsubfig(ivar),AP);   
%         end
    end
end
return