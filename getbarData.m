hfig = gcf;
hsubfig = findobj(hfig,'type','axes');
if ~isempty(hsubfig)
    for is = 1:length(hsubfig)
        hbar = findobj(hsubfig(is),'type','bar');
        Y = get(hbar,'YData');
    end
end
        
        XD = zeros(nmodels,length(hbar));
        YD = zeros(nmodels,length(hbar));
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
for j = 1:size(XD,2)
    figure
    plot(XD(:,j),YD(:,j),'LineStyle','none','Marker','o')
%     setPropPlot
end