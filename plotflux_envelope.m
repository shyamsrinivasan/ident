function plotflux_envelope(model,flux,fluxind,hsubfig,prxnid,call)
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
if nargin < 6 || call == 1
    ColorSpec = chooseColors(nflux,{'Red'});
elseif call == 2
    ColorSpec = chooseColors(nflux,{'LightGray'});
end

% hline = zeros(nflux,1);
%Scatter Plot all Fluxes
% fig_name = texlabel(['Flux Envelope \mu = ' num2str(model.gmax) 'h^{-1}']);
fig_name = sprintf('Flux Envelope %g',model.gmax);
for ivar = 1:nflux
    if isfield(model,'Regulators');
        tfp = find(strcmpi(fluxind{ivar},model.Regulators));
        y_label1 = sprintf('Flux %s',model.Regulators{tfp});
    else
        tfp = find(strcmpi(fluxind{ivar},model.rxns));
        y_label1 = sprintf('Flux %s',model.rxns{tfp});
    end
    if any(tfp)
        if isempty(findobj('type','figure','Name',fig_name))
            hfig = figure('Name',fig_name); 
            figure(hfig);
        else
            hfig = findobj('type','figure','Name',fig_name);
            figure(hfig);
        end
        if tfp ~= model.Vupind
            if hsubfig(tfp) ~= 0
                hca = findobj(hsubfig(tfp),'type','axes');
                set(hfig,'CurrentAxes',hca);  
            else%subplot is unassigned
                hsubfig(tfp) = subplot(2,n,ivar);       
                hca = gca;
                %Make sure more plots can be added at end of loop
                set(hca,'NextPlot','add');     
                set(hsubfig(tfp),'NextPlot','add');
            end
%             hsubfig(tfp-1) = subplot(nflux,1,ivar);
            fluxd = fluxt(:,tfp);%flux of ivar
            fluxd_p = fluxt(:,prxnid);
            Y = zeros(length(fluxd),1);
            initf = fluxd(1);
            fluxd(1) = [];
            Y(1) = fluxd_p(1);
            Y(2:length(fluxd)+1) = fluxd_p(2:length(fluxd)+1);
            hline = line([initf;fluxd],Y);
            set(hline,'LineStyle','none',...
                         'Marker','o',...
                         'MarkerSize',5,...
                         'MarkerFaceColor',ColorSpec{ivar},...
                         'MarkerEdgeColor',ColorSpec{ivar});

%             set(get(gca,'YLabel'),'String',y_label1);  
            set(get(gca,'YLabel'),'FontName','CMU Serif');   
            set(get(gca,'YLabel'),'FontSize',10);  

    %         xlabel = sprintf('Steady States');
    %         set(get(gca,'XLabel'),'String',xlabel);  
    %         set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
    %         set(get(gca,'XLabel'),'FontSize',12);
        end
    end
end

for ifig = 1:length(hsubfig)
    if hsubfig(ifig) ~= 0
        hline = findobj(hsubfig(ifig),'type','line');
        nlines = length(hline);
        maxX = zeros(nlines,1);
        minX = zeros(nlines,1);
        maxY = zeros(nlines,1);
        minY = zeros(nlines,1);
        for iline = 1:nlines
            XData = get(hline(iline),'XData');
            YData = get(hline(iline),'YData');
            maxX(iline) = max(XData);
            minX(iline) = min(XData);
            maxY(iline) = max(YData);
            minY(iline) = min(YData);
        end
        maxX = max(maxX);
        minX = min(minX);
        maxY = max(maxY);
        minY = min(minY);
        [minY,maxY] = fixaxisbounds(minY,maxY);
        [minX,maxX] = fixaxisbounds(minX,maxX);
        set(hfig,'CurrentAxes',hsubfig(ifig));
        set(gca,'YLim',[minY maxY]);
        AP.XColor = [.1 .1 .1];
        AP.YColor = [.1 .1 .1];
        AP.XLim = [minX maxX];
%         AP.XTick = 0:1:max(X);
        set(hsubfig(ifig),AP); 
    end
end


return