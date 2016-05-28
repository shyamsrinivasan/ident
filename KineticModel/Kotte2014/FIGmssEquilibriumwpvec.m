function colorSpec =...
FIGmssEquilibriumwpvec(yeq,pvec,allnss,idx,idp,datatype,plotype,colorSpec)
% plot the equilibrium state variable change wrt chage in parameters for
% systems whose multistationarity has been evaluated using MATCONT
% call function with relevant inputs (flux/concentration) to plot figures
if nargin<7
    plotype = 2;
end
if nargin < 8
    colorSpec = [];
end

switch datatype
    case 'flux'
        datatype = 1;
    case 'conc'
        datatype = 2;
    otherwise
        datatype = 0;
end

if plotype == 2
    ndp = length(idp);
elseif plotype == 3
    ndp = 1;
end    

nss = zeros(ndp,1);
for ip = 1:ndp
    % identify number of steady states in y for each parameter set       
    nss(ip) = length(unique(allnss.(['iid' num2str(ip)])));   
end

% maximum number of mss
nss = max(nss);
if isempty(colorSpec)||size(colorSpec,1)<nss
    colorSpec = chooseColors(nss);
end

switch plotype
    case 2
        % 2-D plots
        plot2D(yeq,pvec,allnss,nss,colorSpec,idx,idp,datatype)
    case 3        
        % 3-D plots
        plot3D(yeq,pvec,allnss,nss,colorSpec,idx,idp,datatype)
    otherwise
end
end

function plot3D(yeq,pvec,allnss,nss,colorSpec,idx,idp,datatype)
% 3-D plots : 1 decision variable vs 2 parameters

ndx = length(idx);
ssid = allnss.iid1;
mss = unique(ssid);
% several combinations of 2 parameters
ncomb = nchoosek(idp,2);
for icomb = 1:size(ncomb,1);
    for iy = 1:ndx
        figure
        ydata = yeq(idx(iy),:,1);
        pdata = pvec(:,ncomb(icomb,:),1);
        hlss = zeros(nss,1);
        Point = struct([]);
        for iss = 1:nss
            % line properties
            Point(iss).MarkerEdgeColor = colorSpec{iss};
            Point(iss).MarkerFaceColor = [1 1 1]; % colorSpec{iss};
            Point(iss).MarkerSize = 10;
            Point(iss).Marker = '.';
            Point(iss).LineStyle = 'none'; 
            % plot data as per # mss
            hlss(iss) = plot3(pdata(ssid==mss(iss),1),...
                             pdata(ssid==mss(iss),2),...
                             ydata(ssid==mss(iss))); 
            set(gca,'NextPlot','add');
            if mss(iss)~=0
                ptname = ['Possible LP = ' num2str(mss(iss)-2)];
            else
                ptname = ['Possible LP = ' num2str(0)];
            end
            Point(iss).DisplayName = ptname;
            set(hlss(iss),Point(iss));
        end
        % axes labels
        [xlabel,ylabel,zlabel] =...
        getaxislabels(3,datatype,idx(iy),ncomb(icomb,:));
    
        % axes protperties
        setproperties(3,gcf,xlabel,ylabel,zlabel);
    end
end
end

function plot2D(yeq,pvec,allnss,nss,colorSpec,idx,idp,datatype)
% 2-D plots 1 decision variable vs 1 parameter

figure
ndx = length(idx);
ndp = length(idp);

nfig = ndx*ndp;
if rem(nfig,2)~=0
    nfig = nfig+1;
end
if nfig ~=2
    nrows = nfig/2;
else
    nrows = 2;
end
ncol = 2;

ylabel = cell(nfig,1);
xlabel = cell(nfig,1);
ifig = 1;
for ip = 1:ndp    
    ssid = allnss.(['iid' num2str(ip)]);
    mss = unique(allnss.(['iid' num2str(ip)]));
    for iy = 1:ndx
        % gather data
        ydata = yeq(idx(iy),:,ip);
        pdata = pvec(:,idp(ip),ip);
        hsfig = subplot(nrows,ncol,ifig);
        set(hsfig,'NextPlot','add');
        hlss = zeros(nss,1);
        Point = struct([]);
        for iss = 1:nss
            % line properties
            Point(iss).MarkerEdgeColor = colorSpec{iss};
            Point(iss).MarkerFaceColor = [1 1 1]; % colorSpec{iss};
            Point(iss).MarkerSize = 10;
            Point(iss).Marker = 'o';
            Point(iss).LineStyle = 'none'; 
            % plot data as per # mss
            hlss(iss) = line(pdata(ssid==mss(iss)),ydata(ssid==mss(iss))); 
            if mss(iss)~=0
                ptname = ['Possible LP = ' num2str(mss(iss)-2)];
            else
                ptname = ['Possible LP = ' num2str(0)];
            end
            Point(iss).DisplayName = ptname;
            set(hlss(iss),Point(iss));
        end
        [xlabel{ifig},ylabel{ifig}] =...
        getaxislabels(2,datatype,idx(iy),idp(ip));
        
        % axes protperties
        setproperties(2,hsfig,xlabel{ifig},ylabel{ifig});
        ifig = ifig+1;
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
set(get(gca,'YLabel'),'FontSize',18); 
set(get(gca,'XLabel'),'String',xlabel);  
set(get(gca,'XLabel'),'FontName','Arial');   
set(get(gca,'XLabel'),'FontSize',18);
if plotype == 3
    set(get(gca,'ZLabel'),'String',zlabel);  
    set(get(gca,'ZLabel'),'FontName','Arial');   
    set(get(gca,'ZLabel'),'FontSize',18);
    axesP.ZColor = [.1 .1 .1];
end
axesP.FontName  = 'Arial';
axesP.FontSize = 18;
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

function [xlabel,ylabel,zlabel] = getaxislabels(plotype,datatype,idx,idp)

fluxlist = {'ACpts mmole/h','ENZC mmole/h','ECbiomass(FDP) mmole/h',...
           'GLUX mmole/h','PEPout mmole/h'};
cnclist = {'PEP','FDP','ENZ'};  
parlist = {'Parameter 1, kEcat','','','Parameter 2, vFbpmax','','',...
           'Parameter 3, vEXmax','','','','','','','Parameter 4, kPEPout'};

if plotype == 2
    if datatype == 1
        ylabel = fluxlist(idx);
    elseif datatype == 2
        ylabel = cnclist(idx);
    end
elseif plotype == 3
    ylabel = parlist(idp(2));
    if datatype == 1
        zlabel = fluxlist(idx);
    elseif datatype == 2
        zlabel = cnclist(idx);
    end
end
xlabel = parlist(idp(1));
end

