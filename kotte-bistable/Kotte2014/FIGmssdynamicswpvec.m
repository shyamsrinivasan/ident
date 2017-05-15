function [hfig,hsfig,hline] =...
FIGmssdynamicswpvec(y,tout,pvec,idx,idp,jpts,datatype,tspan,Line,hfig,hsfig,hline)
% function to plot dynamic time profiles of flux/concentration for systems
% whose multistationarity has been evaluated using MATCONT
if nargin <12
    hline = [];
end
if nargin<11
    hsfig = [];
end
if nargin<10
    hfig = [];
end
if nargin<9
    Line = struct([]);
end
if nargin < 8 || isempty(tspan)
    tspan = 1:length(tout);
end
if nargin < 7
    datatype = 2;
end
if nargin < 6
    jpts = 1;
end

if ischar(datatype)
    switch datatype
        case 'flux'
            datatype = 1;
        case 'conc'
            datatype = 2;
        otherwise
            datatype = 0;
    end
end
    

% distinguish between multiple sampled datasets and individual datasets
% suplied in y
if length(size(y))>2
    % goto multidimensional plot as earlier
    plot4Dwithpoints(y,tout,pvec,idx,tspan,idp,jpts,datatype,Line)
else
    % go to single idemnsional plots
    [hfig,hsfig,hline] =...
    plot2Dwopoints(hfig,hsfig,hline,y,tout,idx,tspan,Line,datatype);
end
end

function [] =...
plot4Dwithpoints(y,tout,allpvec,idx,tspan,idp,jpts,datatype,Line)
ndp = length(idp);
ndx = length(idx);

nfig = ndx*ndp;
if rem(nfig,2)~=0
    nfig = nfig+1;
end
if nfig ~=2
    nrows = nfig/2;
else
    nrows = 2;
end
if ndp > 1
    ncol = 2;
else
    ncol=1;
end

figure
ifig = 1;
for ip = 1:ndp
    % sort data as per increasing pvec
    pvec = allpvec(:,:,ip);
    [~,index] = sortrows(pvec,idp(ip));
    id = ismember(index,jpts);
    select_sorted_p = pvec(index(id),:);
    pval = select_sorted_p(:,idp(ip));
    for iy = 1:ndx
        ydata =...
        reshape(y(idx(iy),tspan,index(id),ip),length(tspan),length(jpts));
        pvec_select = pvec(index(id),:);
        xdata = repmat(tout(tspan),1,length(jpts));
        hsfig = subplot(nrows,ncol,ifig);
        set(hsfig,'NextPlot','add');
        ptname = cell(length(jpts),1);
        for ipts = 1:length(jpts)   
            ptname{ipts} = [num2str(jpts(ipts)) '=' num2str(pval(ipts))];             
        end           
        ht = line(xdata,ydata);
        text(xdata(end,:),ydata(end,:),ptname);
        set(ht,{'DisplayName'},ptname);
        
        [xlabel,ylabel] = getaxislabels(2,datatype,idx(iy));              
        setproperties(2,hsfig,xlabel,ylabel);        
        ifig = ifig+1;
    end    
end
end

function [xlabel,ylabel] = getaxislabels(plotype,datatype,idx)

fluxlist = {'ACpts mmole/h','ENZC mmole/h','ECbiomass(FDP) mmole/h',...
           'GLUX mmole/h','PEPout mmole/h'};
cnclist = {'PEP','FDP','ENZ'};  
% parlist = {'Parameter 1, kEcat','','','Parameter 2, vFbpmax','','',...
%            'Parameter 3, vEXmax','','','','','','','Parameter 4, kPEPout'};

if plotype == 2
    if datatype == 1
        ylabel = fluxlist(idx);
    elseif datatype == 2
        ylabel = cnclist(idx);
    end
% elseif plotype == 3
%     ylabel = parlist(idp(2));
%     if datatype == 1
%         zlabel = fluxlist(idx);
%     elseif datatype == 2
%         zlabel = cnclist(idx);
%     end
end
xlabel = sprintf('Time (s)'); 
% xlabel = parlist(idp(1));
end

function setproperties(plotype,hsfig,xlabel,ylabel,zlabel)
if nargin < 5
    zlabel = {};
end
% set current axes
% axes(hsfig);
% set axis properties
set(get(hsfig,'YLabel'),'String',ylabel);  
set(get(hsfig,'YLabel'),'FontName','Arial');   
set(get(hsfig,'YLabel'),'FontSize',22); 
set(get(hsfig,'XLabel'),'String',xlabel);  
set(get(hsfig,'XLabel'),'FontName','Arial');   
set(get(hsfig,'XLabel'),'FontSize',22);
if plotype == 3
    set(get(hsfig,'ZLabel'),'String',zlabel);  
    set(get(hsfig,'ZLabel'),'FontName','Arial');   
    set(get(hsfig,'ZLabel'),'FontSize',22);
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

function [hfig,hsfig,hln] =...
plot2Dwopoints(hfig,hsfig,hln,y,tout,idx,tspan,Line,datatype)
if isempty(hfig)
    hfig = figure;
end
nfig = length(idx);
if rem(nfig,2)~=0
    nfig = nfig+1;
end
% if nfig ~=2
%     nrows = nfig/2;
%     ncol = 2;
% else
%     nrows = 2;
%     ncol = 1;
% end
nrows = nfig;
ncol = 1;

ifig = 1;
if isempty(hsfig)
    hsfig = zeros(nfig,1);
end
for iy = 1:length(idx)
    ydata = y(idx(iy),tspan);
    xdata = tout(tspan);
    if hsfig(ifig) == 0
        hsfig(ifig) = subplot(nrows,ncol,ifig);
    else
%         axes(hsfig(ifig));
        set(hfig,'CurrentAxes',hsfig(ifig));
    end
    set(hsfig(ifig),'NextPlot','add');
%     Line.Color = [0 0 0];
%     Line.DisplayName
    Line.LineWidth = 2.0;
%     Line.LineStyle = ':';
    
    %plot data
    hlss = line(xdata,ydata);
    set(hlss,Line);
    hln = [hln;hlss];   
    
    ifig = ifig+1;
end

for iy = 1:length(idx)
    [xlabel,ylabel] = getaxislabels(2,datatype,idx(iy));
    setproperties(2,hsfig(iy),xlabel,ylabel)
end
end



