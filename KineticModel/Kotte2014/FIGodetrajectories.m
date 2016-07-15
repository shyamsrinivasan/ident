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
    
    [xlabel,ylabel] = getKotteaxislabels(2,datatype,idx);
    setKotteproperties(2,ha,xlabel,ylabel)
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

    [xlabel,ylabel,zlabel] = getKotteaxislabels(3,datatype,idx);
    setKotteproperties(3,ha,xlabel,ylabel,zlabel)
end
end




