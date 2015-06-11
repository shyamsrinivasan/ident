%Plot flux statistical distribution
function plotfluxdistr(model,flux,fluxind,variable)
if nargin < 3
    nflux = length(model.Vind)+length(model.Vexind);
else
    nflux = length(fluxind);
end
hdsubfig = zeros(nflux,1);
if rem(nflux,2)==0
    n = nflux/2;
else
    n = (nflux+1)/2;
end
%flux %row - fluxes %column - samples
fluxt = flux';%row - samples %column - fluxes
maxflux = max(fluxt);
minflux = min(fluxt);
hdbar = zeros(nflux,1);
fColor = [];

switch variable
    case 'concentration'
        Var = model.Metabolites;
    case 'flux'
        if isfield(model,'Regulators')
            Var = model.Regulators;
        else
            Var = model.Enzyme;
        end
end

for i = 1:nflux
    tf = strcmpi(fluxind{i},Var);
    switch variable
        case 'concentration'
            if any(tf)
                y_label1 = sprintf('%s mmole/gDCW',fluxind{i});
                if isempty(findobj('type','figure',...
                                   'Name','Concentration Distribution'))
                    hdfig = figure('Name','Concentration Distribution'); 
                else
                    hdfig = findobj('type','figure','Name','Concentration Distribution');
                end
                figure(hdfig);
                hdsubfig =...
                barplot(i,n,hdfig,hdsubfig,fluxt,tf,y_label1,fColor);
            end
        case 'flux'
            if any(tf)
                y_label1 = sprintf('%s Flux mmole/gDCW/s',fluxind{i});
                if isempty(findobj('type','figure',...
                                   'Name','Flux Distribution'))
                    hdfig = figure('Name','Flux Distribution'); 
                else
                    hdfig = findobj('type','figure','Name','Flux Distribution');
                end
                figure(hdfig);
                hdsubfig =...
                barplot(i,n,hdfig,hdsubfig,fluxt,tf,y_label1,fColor);
            end
    end
end
% for i = 1:nflux
%     if isfield(model,'Regulators')
%         tfp = find(strcmpi(fluxind{i},model.Regulators));
%     elseif isfield(model,'Enzyme')
%         tfp = find(strcmpi(fluxind{i},model.Enzyme));
%     end
%     
%     if any(tfp)
%         if isempty(findobj('type','figure','Name','Flux Distribution'))
%             hdfig = figure('Name','Flux Distribution'); 
%         else
%             hdfig = findobj('type','figure','Name','Flux Distribution');
%         end
%         figure(hdfig);
%         hdsubfig(i) = subplot(2,n,i);  
%         if isempty(findobj(hdsubfig(i),'type','axes'))
%             haxis = axes;
%         else
%             haxis = findobj(hdsubfig(i),'type','axes');
%             set(hdfig,'CurrentAxes',haxis);
%         end
% %         [nelem,brange] = hist(fluxt(:,tfp),2);
%         binEdge = linspace(minflux(tfp),maxflux(tfp),2+1);   
%         [bincounts,ind]= histc(fluxt(:,tfp),[binEdge(1:end-1) Inf]);
%         %Assign Scatter X and Y
%         [X,Y,Xlbl,fColor] = scatterXY(fluxt(:,tfp),bincounts,ind,fColor);
%         hline = bar(X,Y,'BarWidth',0.4,...
%                         'EdgeColor',fColor(1,:),...
%                         'FaceColor',fColor(1,:));       
%         set(gca,'XTickLabel',Xlbl);
% %         hline = line(X,Y,'MarkerFaceColor',fColor(i,:),...
% %                          'MarkerEdgeColor',fColor(i,:));
% %         set(hline,'LineStyle','none',...
% %                   'Marker','o',...
% %                   'MarkerSize',5);            
% %         preass = 0;
% %         hdbar(i) = scatter(X,Y,40,fColor,...
% %                                   'LineWidth',1.5);
% %         set(hdbar(i),'MarkerFaceColor',fColor);
% %         hdbar(i) = bar(haxis,brange,nelem);
%         AP.XColor = [.1 .1 .1];
%         AP.YColor = [.1 .1 .1];
% %         AP.XLim = [0 max(X)+1];
% %         AP.XTick = 0:1:max(X);
% %         miny = min(Y);
% %         maxy = max(Y);
% %         if maxy-miny <= 1e-5
% %             AP.YLim = [0.9*miny 1.1*maxy];
% %         else
% %             AP.YLim = [miny maxy];
% %         end
%         set(hdsubfig(i),AP);
%         ylabel = sprintf('%s Flux mmole/gDCW/s',fluxind{i});
%         set(get(gca,'YLabel'),'String',ylabel);  
%         set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
%         set(get(gca,'YLabel'),'FontSize',12);  
%         xlabel = sprintf('Steady States');
%         set(get(gca,'XLabel'),'String',xlabel);  
%         set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
%         set(get(gca,'XLabel'),'FontSize',12);
% %         set(hdbar(i),'FaceColor',ColorSpec{i});
%     else
%         fprintf('No flux for %s exists',fluxind{i});
%         continue        
%     end
% end

function [X,Y,Xlbl,fColor] = scatterXY(flux,bincounts,ind,fColor)
nsmpl = length(bincounts);
%Randomly choose colors
assign = 1;
ColorSpec = chooseColors(nsmpl);

allind = unique(ind);
% nstates = length(allind);       
% Y = zeros(sum(bincounts),1);
% X = zeros(sum(bincounts),1);
if isempty(fColor)
%     fColor = cell(sum(bincounts),1);
    fColor = cell(length(bincounts),1);
end
% idx = 0;
% for i = 1:length(bincounts)
% % end
% % for i = 1:nstates
%     Y(idx+1:idx+bincounts(i)) = flux(ind == i);
%     X(idx+1:idx+bincounts(i)) = i;
%     if assign
%         fColor(idx+1:idx+bincounts(i)) = ColorSpec(idx+1:idx+bincounts(i));
%     else
%         fColor(idx+1:idx+bincounts(i),:) = ColorSpec(idx+1:idx+bincounts(i),:);
%     end
%     idx = idx + bincounts(i);        
% end
X = zeros(length(bincounts),1);
Y = zeros(length(bincounts),1);
Xlbl = zeros(length(bincounts),1);
if length(allind) > 1
    for i = 1:length(bincounts)
        X(i) = i;
        try
            fColor(i,:) = ColorSpec{i};
        catch
            fColor(i,:) = ColorSpec(i,:);
        end
        if ~isempty(flux(ind==i))
            if length(unique(flux(ind==i))) == 1

                Xlbl(i) = unique(flux(ind == i));
            else
                Xs = unique(flux(ind==i));
                Xlbl(i) = Xs(1);
            end
            Y(i) = bincounts(i);
        else
            X(i) = 0;
            Y(i) = 0;
        end    
    end
else
    X(1) = 1;
    try
        fColor(1,:) = ColorSpec{1};
    catch
        fColor(1,:) = ColorSpec(1,:);
    end
    Y(1) = sum(bincounts);
    Xlbl(1) = flux(1);
end
X(Y==0) = [];
Xlbl(Y==0)=[];
fColor(Y==0,:) = [];
Y(Y==0) = [];
Xlbl = cellstr(num2str(Xlbl));
try
    fColor = cell2mat(fColor);
catch
    fColor = fColor;
end

function hsubfig = barplot(j,nplots,hfig,hsubfig,valt,tf,y_label1,fColor)
    maxval = max(valt);
    minval = min(valt);
    if hsubfig(j) ~= 0
        hca = findobj(hsubfig(j),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(j) = subplot(2,nplots,j);       
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(j),'NextPlot','add');
    end
    binEdge = linspace(minval(tf),maxval(tf),2+1);
    [bincounts,ind]= histc(valt(:,tf),[binEdge(1:end-1) Inf]);
    
    %Assign Scatter X and Y
    [X,Y,Xlbl,fColor] = scatterXY(valt(:,tf),bincounts,ind,fColor);
    hbar = bar(X,Y,'BarWidth',0.4,...
                    'EdgeColor',fColor(1,:),...
                    'FaceColor',fColor(1,:));       
    set(gca,'XTickLabel',Xlbl);
    
    AP.XColor = [.1 .1 .1];
    AP.YColor = [.1 .1 .1];
    
    set(hsubfig(j),AP);
%     ylabel = sprintf('%s Flux mmole/gDCW/s',fluxind{i});
    set(get(gca,'YLabel'),'String',y_label1);  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'YLabel'),'FontSize',12);  
    xlabel = sprintf('Steady States');
    set(get(gca,'XLabel'),'String',xlabel);  
    set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'XLabel'),'FontSize',12);
return
