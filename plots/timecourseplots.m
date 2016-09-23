function [hfig,ha] = timecourseplots(xdata,ally,datatype,name,lcmodel,hfig,ha,LP,addanot)
if nargin<9
    addanot = [];
end
if nargin<8
    LP = [];
end
if nargin<7
    ha = [];
end
if nargin<6
    hfig = [];
end

global Line annot type model
Line = LP;
annot = addanot;
type = datatype;
model = lcmodel;
if isempty(Line)
    Line.Color = [0 0 0];
    Line.LineWidth = 2;
end

switch(type)
    case 1        
        % get id for metabolites in cell array name
        if iscell(name)
            name = cellfun(@(x)strcmpi(model.mets,x),name,'UniformOutput',false);
            name = cellfun(@(x)find(x),name,'UniformOutput',false);
            name = name(~cellfun('isempty',name));
            yid = cell2mat(name);
        else
            yid = name;
        end    
        % get ydata and plot info
        ydata = ally(yid,:);
        [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,yid,hfig,ha);
        drawnow
    case 2
        % get id for fluxes in cell array name
        if iscell(name)
            name = cellfun(@(x)strcmpi(model.rxns,x),name,'UniformOutput',false);
            name = cellfun(@(x)find(x),name,'UniformOutput',false);
            name = name(~cellfun('isempty',name));
            yid = cell2mat(name);
        else
            yid = name;
        end  
        % get ydata and plot info
        ydata = ally(yid,:);
        [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,yid,hfig,ha); 
        drawnow
    otherwise
        fprintf('No suitable plot functions found\n');
end
end

function [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,idx,hfig,ha)
global Line annot type model
if ~isempty(hfig)
    set(0,'CurrentFigure',hfig);
else
    hfig = figure('Name','Time Courses'); 
    set(0,'CurrentFigure',hfig);
end

ndata = size(ydata,1);
if ndata>1    
    if rem(ndata,2)==0
        nrows = ndata/2;
    else
        nrows = (ndata+1)/2;
    end    
else
    nrows = 1;
end
ncols = 2;
nplots = nrows*ncols;
if isempty(ha)|length(ha)~=nplots
    ha = zeros(nplots,1);
end

for ip = 1:nplots
    if ip<=ndata
        if ha(ip)   % if subplot is already assigned
            set(hfig,'CurrentAxes',ha(ip));
            hlval = line(xdata,ydata(ip,:));
        else        % if subplot is emtpy and unassigned
            ha(ip) = subplot(nrows,ncols,ip);
            set(hfig,'CurrentAxes',ha(ip));
            hlval = line(xdata,ydata(ip,:));
        end
        if ~isempty(Line)
            set(hlval,Line);
        end
        % set axis labels
        if type==1
            [xlabel,ylabel] = getaxislabels(2,11,model,[1 idx(ip)]);
        elseif type==2
            [xlabel,ylabel] = getaxislabels(2,12,model,[1 idx(ip)]);
        end
        % set axis properties
        setproperties(2,ha(ip),xlabel,ylabel);
    else
        break
    end
end
end
