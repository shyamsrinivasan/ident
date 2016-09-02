function [hfig,ha] = timecourseplots(xdata,ally,type,name,model,hfig,ha,LP,addanot)
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

global Line annot
Line = LP;
annot = addanot;

switch(type)
    case 1        
        % get id for metabolites in cell array name
        if iscell(name)
            name = cellfun(@(x)strcmpi(model.mets,x),name,'UniformOutput',false);
            name = cellfun(@(x)find(x),name,'UniformOutput',false);
            yid = cell2mat(name);
        else
            yid = name;
        end    
        % get ydata and plot info
        ydata = ally(yid,:);
        [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,hfig,ha,LP);
        drawnow
    case 2
        % get id for fluxes in cell array name
        if iscell(name)
            name = cellfun(@(x)strcmpi(model.rxns,x),name,'UniformOutput',false);
            name = cellfun(@(x)find(x),name,'UniformOutput',false);
            yid = cell2mat(name);
        else
            yid = name;
        end  
        % get ydata and plot info
        ydata = ally(yid,:);
        [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,hfig,ha,LP); 
        drawnow
    otherwise
        fprintf('No suitable plot functions found\n');
end
end

function [hfig,ha,hlval] = plot2Dcourse(xdata,ydata,hfig,ha)
global Line annot
if ~isempty(hfig)
    set(0,'CurrentFigure',hfig);
else
    hfig = figure('Name','Met Courses'); 
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
        % set axis properties
    else
        break
    end
end
end
