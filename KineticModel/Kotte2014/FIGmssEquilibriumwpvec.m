function colorSpec =...
FIGmssEquilibriumwpvec(yeq,pvec,allnss,idx,idp,datatype,colorSpec)
% plot the equilibrium state variable change wrt chage in parameters for
% systems whose multistationarity has been evaluated using MATCONT
% call function with relevant inputs (flux/concentration) to plot figures
if nargin < 7
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

figure
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
ncol = 2;
% if ndp > 1
%     ncol = 2;
% else
%     ncol=1;
% end

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
        if datatype == 1
            switch idx(iy)
                case 1
                    ylabel{ifig} = sprintf('ACpts mmole/h');
                case 2
                    ylabel{ifig} = sprintf('ENZC mmole/h');
                case 3
                    ylabel{ifig} = sprintf('ECbiomass(FDP) mmole/h');
                case 4
                    ylabel{ifig} = sprintf('GLUX mmole/h');
                case 5
                    ylabel{ifig} = sprintf('PEPout mmole/h');
            end    
        elseif datatype == 2
            switch idx(iy)
                case 1
                    ylabel{ifig} = sprintf('PEP');
                case 2
                    ylabel{ifig} = sprintf('FDP');
                case 3
                    ylabel{ifig} = sprintf('ENZ');                
            end 
        end
        switch idp(ip)
            case 1
                xlabel{ifig} = 'Parameter 1, kEcat';
            case 4
                xlabel{ifig} = 'Parameter 2, vFbpmax';
            case 7
                xlabel{ifig} = 'Parameter 3, vEXmax';
            case 14
                xlabel{ifig} = 'Parameter 4, kPEPout';
        end
        % axes protperties
        set(get(gca,'YLabel'),'String',ylabel{ifig});  
        set(get(gca,'YLabel'),'FontName','Arial');   
        set(get(gca,'YLabel'),'FontSize',18); 
        set(get(gca,'XLabel'),'String',xlabel{ifig});  
        set(get(gca,'XLabel'),'FontName','Arial');   
        set(get(gca,'XLabel'),'FontSize',18);
        axesP.FontName  = 'Arial';
        axesP.FontSize = 18;
        axesP.LineWidth = 1.5;
        axesP.TickLength = [0.01 0.01];
        axesP.XColor = [.1 .1 .1];
        axesP.YColor = [.1 .1 .1];
        set(hsfig,axesP);
        ifig = ifig+1;
    end    
end

% set figure protperties, axes, etc

    
