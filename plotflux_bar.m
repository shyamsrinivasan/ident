function plotflux_bar(model,flux,fluxind)
if nargin < 3
    nflux = length(model.Vind)+length(model.Vexind);
else
    nflux = length(fluxind);
end
hsubfig = zeros(nflux,1);
nsmpl = size(flux,2);

%Randomly choose colors
ColorSpec = chooseColors(nflux);

%Bins for all fluxes
%flux %row - fluxes %column - samples
fluxt = flux';%row - samples %column - fluxes
fluxiqr = iqr(fluxt);
maxflux = max(fluxt);
minflux = min(fluxt);
h = 2.*fluxiqr.*(nsmpl^(-1/3));
nbins = (maxflux-minflux)./h;
hbar = zeros(nflux,1);
%Flux Plot (Direct plot of fluxes on bar graph for all samples)
for i = 1:nflux
    if isfield(model,'Regulators');
        tfp = find(strcmpi(fluxind{i},model.Regulators));
    else
        tfp = find(strcmpi(fluxind{i},model.rxns));
    end
    if any(tfp)
%         xticklabel{i} = sprintf('%s Flux',fluxind{i});
        if isempty(findobj('type','figure','Name','Fluxes'))
            hfig = figure('Name','Fluxes'); 
            figure(hfig);
        end 
        figure(hfig);
        hsubfig(i) = subplot(nflux,1,i);
        fluxd = flux(tfp,:);
        hbar(i) = bar(hsubfig(i),1:length(fluxd),fluxd,'BarWidth',0.4,...
                                                       'EdgeColor',ColorSpec{i},...
                                                       'FaceColor',ColorSpec{i});
        hca = findobj(hsubfig(i),'type','axes');
        set(hfig,'CurrentAxes',hca);
%         set(gca,'Title',text('String',pfname,'Color','k'));
        set(get(gca,'Title'),'FontName','CMU Serif');
        set(get(gca,'Title'),'FontSize',18);
        set(get(gca,'Title'),'FontWeight','bold'); 
        ylabel = sprintf('%s Flux',fluxind{i});
        set(get(gca,'YLabel'),'String',ylabel);  
        set(get(gca,'YLabel'),'FontName','CMU Serif');   
        set(get(gca,'YLabel'),'FontSize',18);    
%         set(hsubfig(i),'Ylim',[min(min(fluxd)),max(max(fluxd))]);
        xticklabel = [{'Initial'} num2cell(1:length(fluxd)-1)];
%         set(hsubfig(i),'XTickLabel','none');
        set(hbar(i),'FaceColor',ColorSpec{i});
%         set(hbar(i),'EdgeColor',ColorSpec{i});
%         for ibar = 1:nsmpl
%             set(hbar(i,ibar),'FaceColor',ColorSpec{ibar});
%         end
        hold on
    else
        fprintf('No flux for %s exists',fluxind{i});
        continue        
    end   
end

return