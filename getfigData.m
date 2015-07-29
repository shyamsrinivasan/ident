hfig = gcf;
hsubfig = findobj(hfig,'type','line');
XD = zeros(nmodels,length(hsubfig));
YD = zeros(nmodels,length(hsubfig));
for i = 1:length(hsubfig)
    X = get(hsubfig(i),'XData');
    if length(X)>200
        XD(:,i) = get(hsubfig(i),'XData');
        YD(:,i) = get(hsubfig(i),'YData');
    end
end
% XD(:,XD==0) = [];
% YD(:,YD==0) = [];
for j = 1:size(YD,2)
    figure
    %plot Flux Distribution wrt other fluxes
%     plot(XD(:,j),YD(:,j),'LineStyle','none','Marker','o')
    %plot concentration distribution wrt other concentrations
    for k = 1:size(YD,2)
        if j ~= k
            plot(YD(:,j),YD(:,k),'LineStyle','none','Marker','o');
        end
    end
%     setPropPlot

end

plot(conc(:,1),conc(:,2),'LineStyle','none','Marker','o');

