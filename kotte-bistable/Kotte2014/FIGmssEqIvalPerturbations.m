function [hfig,ha,hline] =...
FIGmssEqIvalPerturbations(val1,val2,datatype,idx,hfig,ha,Point,addanot)
% val1      - initial value vector
% val2      - equilibrium value vector
% datatype  - 1 for flux, 2 for concentrations
% idx       - index in val1 and val2 to be plotted (row/column vector)
% hfig      - figure handle of existing figure (optional)
% ha        - axis handle of existing figure (optional)
% Point     - structure of point properties (optional)
if nargin<8
    addanot = [];
end
if nargin<7
    Point = [];
end
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
global annot
annot = addanot;

if length(idx) >= 3
    [ha,hline] = plot3DeqPoints(val1,val2,idx,datatype,hfig,ha,Point);
elseif length(idx) >= 2
    [ha,hline] = plot2DeqPoints(val1,val2,idx,datatype,hfig,ha,Point);
end
end

function [ha,hlval] = plot2DeqPoints(val1,val2,idx,datatype,hfig,ha,Point)
% 2D plots of equilibirum points and initial values
global annot
set(0,'CurrentFigure',hfig);
if isempty(ha)
    hlval1 = line(val1(idx(1),:),val1(idx(2),:),...
                'Marker','.',...
                'MarkerSize',8,...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerEdgeColor','k',...
                'LineStyle','none');  
    ha = get(hlval1,'Parent');
    set(ha,'NextPlot','add');
    set(hfig,'CurrentAxes',ha);
    hlval2 = line(val2(idx(1)),val2(idx(2)),...
                'Marker','p',...
                'MarkerSize',16,...            
                'LineStyle','none');    

else
    set(hfig,'CurrentAxes',ha);
    hlval1 = line(val1(idx(1),:),val1(idx(2),:),...
                'Marker','.',...
                'MarkerSize',8,...
                'MarkerFaceColor',[.49 1 .63],...
                'MarkerEdgeColor','k',...
                'LineStyle','none'); 
    hlval2 = line(val2(idx(1)),val2(idx(2)),...
                'Marker','p',...
                'MarkerSize',16,...            
                'LineStyle','none');
end
if ~isempty(Point)
    set(hlval2,Point);
    set(hlval1,Point);
end
% add annotation to data points
if ~isempty(annot)
%     text(val1(idx(1)),val1(idx(2)),annot.text);
    text(val2(idx(1)),val2(idx(2)),annot.text);
end

[xlabel,ylabel] = getKotteaxislabels(2,datatype,idx);

setKotteproperties(3,ha,xlabel,ylabel)
hlval = [hlval1 hlval2];
% set(gac,'XLim',[0 max(xdata)],'YLim',[0 max(ydata)]);
drawnow;
end

function [ha,hlval] = plot3DeqPoints(val1,val2,idx,datatype,hfig,ha,Point)
% 3D plots of equilibrium points and initial values
set(0,'CurrentFigure',hfig);
if isempty(ha)
    hlval1 = line(val1(idx(1),:),val1(idx(2),:),val1(idx(3),:),...
                    'Marker','.',...
                    'MarkerSize',8,...
                    'MarkerFaceColor',[.49 1 .63],...
                    'MarkerEdgeColor','k',...
                    'LineStyle','none');   
%     hlval1 = plot3(val1(idx(1)),val1(idx(2)),val1(idx(3)),...
%                     'Marker','o',...
%                     'MarkerSize',6,...
%                     'MarkerFaceColor',[.49 1 .63],...
%                     'MarkerEdgeColor','k',...
%                     'LineStyle','none');   
    ha = get(hlval1,'Parent');
    set(ha,'NextPlot','add');
    set(hfig,'CurrentAxes',ha);
    hlval2 = line(val2(idx(1)),val2(idx(2)),val2(idx(3)),...
                'Marker','p',...
                'MarkerSize',16,...
                'LineStyle','none');    
    hold on
else
    set(hfig,'CurrentAxes',ha);
    hlval1 = line(val1(idx(1),:),val1(idx(2),:),val1(idx(3),:),...
                    'Marker','.',...
                    'MarkerSize',8,...
                    'MarkerFaceColor',[.49 1 .63],...
                    'MarkerEdgeColor','k',...
                    'LineStyle','none'); 
    hlval2 = line(val2(idx(1)),val2(idx(2)),val2(idx(3)),...
                'Marker','p',...
                'MarkerSize',16,...
                'LineStyle','none');  
%     if ~isempty(Point)
%         hlval2 = plot3(ha,val2(idx(1)),val2(idx(2)),val2(idx(3)),...
%                     'Marker','p',...
%                     'MarkerSize',16,...
%                     'MarkerFaceColor',Point.MarkerEdgeColor,...
%                     'MarkerEdgeColor',Point.MarkerEdgeColor,...
%                     'LineStyle','none');  
%     else
%         
%     end
end
% set point properties
if ~isempty(Point)
    set(hlval2,Point);
    set(hlval1,Point);
end

[xlabel,ylabel,zlabel] = getKotteaxislabels(3,datatype,idx);

setKotteproperties(3,ha,xlabel,ylabel,zlabel)
hlval = [hlval1 hlval2];
drawnow;
end