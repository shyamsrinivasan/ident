%function [hfig,hsubfig] = compareFluxes(model,npertb,varname,solution,...
%                                          ColorSpec,LineP,hsubfig)
%**************************************************************************
%Compare flux response to input perturbations
%December 2014
%**************************************************************************
function [hfig,hsubfig] = compareFluxes(nmodels,npertb,flux_ind,solution,ColorSpec)
nflux = length(flux_ind);
if nargin < 7
    hsubfig = zeros(npertb,1);
end
if nargin < 6
    LineP = struct();
end

flux = zeros(nflux,nmodels);
xticklabel = cell(nflux,1);
DisplayName = cell(nmodels,1);
hsubfig = zeros(npertb,1);
hbar = zeros(npertb,nmodels);
for ipertb = 1:npertb
    pname = sprintf('pertb%d',ipertb);
    pfname = sprintf('Perturbation %d',ipertb);
    for ivar = 1:nflux    
        xticklabel{ivar} = sprintf('Flux %d',flux_ind(ivar));
        for imodel = 1:nmodels
            mname = sprintf('model%d',imodel);
            DisplayName{imodel} = mname;
            flux(ivar,imodel) = solution.(pname).(mname).flux(flux_ind(ivar));
        end       
    end      
    if isempty(findobj('type','figure','Name','Fluxes'))
        hfig = figure('Name','Fluxes'); 
        figure(hfig);
    end  
    figure(hfig);
    hsubfig(ipertb) = subplot(npertb,1,ipertb);     
    hbar(ipertb,:) = bar(hsubfig(ipertb),1:nflux,flux,'BarWidth',1,...
                                            'EdgeColor','none');
    hca = findobj(hsubfig(ipertb),'type','axes');
    set(hfig,'CurrentAxes',hca);
    set(gca,'Title',text('String',pfname,'Color','k'));
    set(get(gca,'Title'),'FontName','Arabic Type Setting');
    set(get(gca,'Title'),'FontSize',12);
    set(get(gca,'Title'),'FontWeight','bold'); 
    set(get(gca,'YLabel'),'String','Fluxes mmole/gDCW.h');  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'YLabel'),'FontSize',12);    
    set(hsubfig(ipertb),'Ylim',[min(min(flux)),max(max(flux))]);
    set(hsubfig(ipertb),'XTickLabel','none');
    for ibar = 1:nmodels
        set(hbar(ipertb,ibar),'FaceColor',ColorSpec{ibar});
    end
    hold on  
end
%Formatting
%Identify Top row
ht_bar = hbar(1,:);
for ih = 1:nmodels
    set(ht_bar(ih),'DisplayName',DisplayName{ih});
end
%Identify middle row
if rem(npertb,2) == 0
    mpertb = npertb/2;
else
    mpertb = (npertb+1)/2;
end
hm_fig = hsubfig(mpertb);
% hca = findobj(hm_bar(1,1),'type','axes');

%Identifying last row
he_fig = hsubfig(end);
% hca = findobj(he_bar(1,1),'type','axes');
set(hfig,'CurrentAxes',he_fig);
set(gca,'XTickLabel',xticklabel);
set(gca,'FontName','Arabic Type Setting');
set(gca,'FontSize',12);
set(get(gca,'XLabel'),'String','Flux #');
set(get(gca,'XLabel'),'FontName','Arabic Type Setting');
set(get(gca,'XLabel'),'FontSize',12);
FProperty.NumberTitle = 'off';
FProperty.Color = [1 1 1];  
set(hfig,FProperty);
legend(findobj(hsubfig(1),'type','axes'),'show');
legend(findobj(hsubfig(1),'type','axes'),'boxoff');

return
    