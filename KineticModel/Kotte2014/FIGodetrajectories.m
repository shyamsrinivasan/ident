function [hfig,ha] =...
FIGodetrajectories(xdyn,ival,xeq,datatype,idx,hfig,ha,Line,Point)
if nargin < 9
    Point = [];
end
if nargin < 8
    Line = [];
end

if length(idx) >= 3
    [hfig,ha] =...
    plot3Dtrajectories(ival,xeq,xdyn,idx,datatype,hfig,ha,Line,Point);
elseif length(idx) >= 2
    plot2Dtrajectories()
end

end

function [heqfig,ha] =...
plot3Dtrajectories(ival,xeq,xdyn,idx,datatype,heqfig,ha,Line,Point)

if length(idx)>3 % get combination of 3 data points
    ncomb = nchoosek(idx,3);
else
    ncomb = idx;
end
for icomb = 1:size(ncomb,1)
    if isempty(heqfig)
        heqfig = figure;
    end
    data = xdyn(ncomb(icomb,:),:);
    xdata = data(1,:);
    ydata = data(2,:);
    zdata = data(3,:);
    
    % plot trajectory
    if isempty(ha)
        hline = plot3(xdata,ydata,zdata);
        ha = get(hline,'Parent');
    else
        hline = plot3(ha,xdata,ydata,zdata);
    end
    hold on
    
    % set line properties
    if ~isempty(Line)
        set(hline,Line);
    end    
    
    % plot the points (initial va;ue and final equilibrium  
    % and set properties
    [heqfig,ha] =...
    FIGmssEqIvalPerturbations(ival,xeq,datatype,ncomb(icomb,:),heqfig,ha,Point);
end
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

function plot2Dtrajectories

end