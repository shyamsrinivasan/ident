function [hfig,ha,hline] =...
FIGodetrajectories(xdyn,ival,xeq,datatype,idx,hfig,ha,line,point)
global Line Point
if nargin < 9
    point = [];
end
if nargin < 8
    line = [];
end
if nargin < 7
    ha = [];
end
if nargin < 6
    hfig = [];
end
Line = line;
Point = point;
if length(idx) >= 3
    [hfig,ha,hline] =...
    plot3Dtrajectories(ival,xeq,xdyn,idx,datatype,hfig,ha);
elseif length(idx) >= 2
    [hfig,ha,hline] =...
    plot2Dtrajectories(ival,xeq,xdyn,idx,datatype,hfig,ha);
end

end

function [heqfig,ha,hline] =...
plot2Dtrajectories(ival,xeq,xdyn,idx,datatype,heqfig,ha)
global Line Point
if length(idx)>2 % get combination of 3 data points
    ncomb = nchoosek(idx,2);
else
    ncomb = idx;
end
for icomb = 1:size(ncomb,1)
    if isempty(heqfig)
        heqfig = figure;
    end
    set(0,'CurrentFigure',heqfig);
    data = xdyn(ncomb(icomb,:),:);
    xdata = data(1,:);
    ydata = data(2,:);
        
    % plot trajectory
    if isempty(ha)
        hline = plot(xdata,ydata);
        ha = get(hline,'Parent');
    else
        hline = plot(ha,xdata,ydata);
    end
    hold on
    
    % set line properties
    if ~isempty(Line)
        set(hline,Line);
    end    
    
    % plot the points (initial va;ue and final equilibrium  
    % and set properties
%     [heqfig,ha] =...
%     FIGmssEqIvalPerturbations(ival,xeq,datatype,ncomb(icomb,:),heqfig,ha,Point);
    
    [xlabel,ylabel] = getaxislabels(2,datatype,idx);
    setproperties(2,ha,xlabel,ylabel)
end

end

function [heqfig,ha,hline] =...
plot3Dtrajectories(ival,xeq,xdyn,idx,datatype,heqfig,ha)
global Line Point
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
%     [heqfig,ha] =...
%     FIGmssEqIvalPerturbations(ival,xeq,datatype,ncomb(icomb,:),heqfig,ha,Point);

    [xlabel,ylabel,zlabel] = getaxislabels(3,datatype,idx);
    setproperties(3,ha,xlabel,ylabel,zlabel)
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
        xlabel = cnclist(idx(1));
        ylabel = cnclist(idx(2));
        zlabel = cnclist(idx(3));
    end
elseif plotype ==2
    if datatype == 1
        [xlabel,ylabel] = deal(fluxlist(idx));
%         zlabel = fluxlist(idx);
%         ylabel = parlist(idp(2));
    elseif datatype == 2
        xlabel = cnclist(idx(1));
        ylabel = cnclist(idx(2));        
    end
    zlabel = {};
end
end

function setproperties(plotype,hsfig,xlabel,ylabel,zlabel)
if nargin < 8
    zdata = [];
end
if nargin < 7
    ydata = [];
end
if nargin < 6
    xdata = [];
end
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
% if ~isempty(xdata)
%     axesP.XLim = [0 max(xdata)];
% end
% if ~isempty(ydata)
%     axesP.YLim = [0 max(ydata)];
% end
% if ~isempty(zdata)
%     axesP.ZLim = [0 max(zdata)];
% end
if plotype == 3
    set(gca,axesP);
else
    set(hsfig,axesP);
end
end

