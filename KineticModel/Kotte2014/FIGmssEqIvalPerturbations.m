function [hfig,ha] =...
FIGmssEqIvalPerturbations(val1,val2,datatype,idx,hfig,ha)
if nargin < 6
    ha = [];
end
if nargin < 5
    hfig = [];
end
if nargin < 4
    idx = [];
end
if nargin < 3
    datatype = 2;
end

if isempty(hfig)
    hfig = figure;
end

if datatype == 2
    if isempty(ha)
        hline = plot3([val1(1) val2(1)],...
                        [val1(2) val2(2)],...
                        [val1(3) val2(3)],...
                        'Marker','o',...
                        'MarkerSize',6,...
                        'MarkerFaceColor',[.49 1 .63],...
                        'MarkerEdgeColor','k');    
        ha = get(hline,'Parent');
        hold on
    else
        plot3(ha,[val1(1) val2(1)],...
                        [val1(2) val2(2)],...
                        [val1(3) val2(3)],...
                        'Marker','o',...
                        'MarkerSize',6,...
                        'MarkerFaceColor',[.49 1 .63],...
                        'MarkerEdgeColor','k'); 
    end
    [xlabel,ylabel,zlabel] = getaxislabels(3,datatype);
elseif datatype == 1
    if isempty(ha)
        hline = plot3([val1(idx(1)) val2(idx(1))],...
                        [val1(idx(2)) val2(idx(2))],...
                        [val1(idx(3)) val2(idx(3))],...
                        'Marker','o',...
                        'MarkerSize',6,...
                        'MarkerFaceColor',[.49 1 .63],...
                        'MarkerEdgeColor','k');    
        ha = get(hline,'Parent');
        hold on
    else
        plot3(ha,[val1(1) val2(1)],...
                        [val1(2) val2(2)],...
                        [val1(3) val2(3)],...
                        'Marker','o',...
                        'MarkerSize',6,...
                        'MarkerFaceColor',[.49 1 .63],...
                        'MarkerEdgeColor','k'); 
    end
    [xlabel,ylabel,zlabel] = getaxislabels(3,datatype,idx);
end   
setproperties(3,ha,xlabel,ylabel,zlabel)

end

function [xlabel,ylabel,zlabel] = getaxislabels(plotype,datatype,idx)
if nargin <3
    idx = [];
end

fluxlist = {'ACpts mmole/h','ENZC mmole/h','ECbiomass(FDP) mmole/h',...
           'GLUX mmole/h','PEPout mmole/h'};
cnclist = {'PEP','FDP','ENZ'};  
% parlist = {'Parameter 1, kEcat','','','Parameter 2, vFbpmax','','',...
%            'Parameter 3, vEXmax','','','','','','','Parameter 4, kPEPout'};

if plotype == 3    
    if datatype == 1
        [xlabel,ylabel,zlabel] = deal(fluxlist(idx));
%         zlabel = fluxlist(idx);
%         ylabel = parlist(idp(2));
    elseif datatype == 2
        xlabel = cnclist(1);
        ylabel = cnclist(2);
        zlabel = cnclist(3);
    end
end
end

function setproperties(plotype,hsfig,xlabel,ylabel,zlabel)
if nargin < 5
    zlabel = {};
end
% set axis properties
set(get(gca,'YLabel'),'String',ylabel);  
set(get(gca,'YLabel'),'FontName','Arial');   
set(get(gca,'YLabel'),'FontSize',22); 
set(get(gca,'XLabel'),'String',xlabel);  
set(get(gca,'XLabel'),'FontName','Arial');   
set(get(gca,'XLabel'),'FontSize',22);
if plotype == 3
    set(get(gca,'ZLabel'),'String',zlabel);  
    set(get(gca,'ZLabel'),'FontName','Arial');   
    set(get(gca,'ZLabel'),'FontSize',22);
    axesP.ZColor = [.1 .1 .1];
end
axesP.FontName  = 'Arial';
axesP.FontSize = 22;
axesP.LineWidth = 1.5;
axesP.TickLength = [0.01 0.01];
axesP.XColor = [.1 .1 .1];
axesP.YColor = [.1 .1 .1];
if plotype == 3
    set(gca,axesP);
else
    set(hsfig,axesP);
end
end