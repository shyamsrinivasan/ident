function bifurcationPlot(y,s1,f1,idx,hfig)
% idx - [x-axis y-axis z-axis] index values
if nargin<7
    hfig = figure;
end

% nvar = size(f1,1);
% npar = size(x1,1)-size(f1,1);

% y = x1(1:nvar,:);
% p = x1(nvar+1:nvar+npar,:);

if size(s1,1)>2
    xindex = cat(1,s1.index);
%     xindex = xindex(2:end-1);
end
evalRe = real(f1);

for kndex = 1:length(xindex)    
    % plot data from initial value to first bifurcation
    if xindex(kndex)+1<=size(y,2)
        if any(evalRe(:,xindex(kndex))<0) && any(evalRe(:,xindex(kndex)+1)>=0)
            % entering unstable from stable region 
            LineP.LineStyle = '--';
            LineP.Color = 'k';
        elseif any(evalRe(:,xindex(kndex))>=0) && any(evalRe(:,xindex(kndex)+1)<0)
            % entering stable from unstable region
            LineP.LineStyle = '-';
            LineP.Color = 'k';
        elseif any(evalRe(:,xindex(kndex))<0) && any(evalRe(:,xindex(kndex)+1)<0)
            % staying in stable region
            LineP.LineStyle = '-';
            LineP.Color = 'k';
        elseif any(evalRe(:,xindex(kndex))>=0) && any(evalRe(:,xindex(kndex)+1)>=0)
            % staying in unstable region
            LineP.LineStyle = '--';
            LineP.Color = 'k';
        end
    end
    LineP.LineWidth = 3;
    LPval = y(:,xindex(2:end-1));
    
    if length(idx)<=2
        % plot2D
        if (kndex+1)<=length(xindex)
            xval = y(idx(1),xindex(kndex):xindex(kndex+1));
            yval = y(idx(2),xindex(kndex):xindex(kndex+1));
        else
            xval = y(idx(1),xindex(kndex):end);
            yval = y(idx(2),xindex(kndex):end);
        end
        % plot data
        hfig = plot2Dbifurcation(xval,yval,idx(1),idx(2),LPval,hfig,LineP);
    elseif length(idx)<=3
        % plot3D 
        if (kndex+1)<=length(xindex)
            xval = y(idx(1),xindex(kndex):xindex(kndex+1));
            yval = y(idx(2),xindex(kndex):xindex(kndex+1));
            zval = y(idx(3),xindex(kndex):xindex(kndex+1));
        else
            xval = y(idx(1),xindex(kndex):end);
            yval = y(idx(2),xindex(kndex):end);
            zval = y(idx(3),xindex(kndex):end);
        end
        % plot data
        hfig = plot3Dbifurcation(xval,yval,zval,idx(1),idx(2),idx(3),...
                                 LPval,hfig,LineP);
    end
end
end

function hfig = plot2Dbifurcation(xval,yval,xid,yid,LPval,hfig,LineP)
if ~isempty(hfig)
    figure(hfig);
    set(gca,'NextPlot','add');
    hl = plot(xval,yval);
else
    hl = plot(xval,yval);
    hfig = gcf;
end

if ~isempty(LineP)
    set(hl,LineP);
else
    set(hl,'LineWidth',3);
end

[xlabel,ylabel] = getaxislabels(2,2,[xid yid]);

% plot bifurcation points
line(LPval(xid,:),LPval(yid,:),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
setproperties(2,gca,xlabel,ylabel);
end

function hfig = plot3Dbifurcation(xval,yval,zval,xid,yid,zid,LPval,hfig,LineP)
if ~isempty(hfig)
    figure(hfig);    
    hl = plot3(xval,yval,zval);
    set(get(hl,'Parent'),'NextPlot','add');
else
    hl = plot3(xval,yval,zval);
    hfig = gcf;
end

if ~isempty(LineP)
    set(hl,LineP);
else
    set(hl,'LineWidth',3);
end

[xlabel,ylabel,zlabel] = getaxislabels(3,2,[xid yid zid]);
plot3(LPval(xid,:),LPval(yid,:),LPval(zid,:),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
setproperties(3,gca,xlabel,ylabel,zlabel);  
end

function [xlabel,ylabel,zlabel] = getaxislabels(plotype,datatype,idx)
if nargin <3
    idx = [];
end

fluxlist = {'ACpts mmole/h','ENZC mmole/h','ECbiomass(FDP) mmole/h',...
           'GLUX mmole/h','PEPout mmole/h'};
cnclist = {'PEP','FDP','ENZ'};  
parlist = {'kEcat','KEacetate','KFbpFBP','vFbpmax','Lfbp','KFbpPEP',...
           'vEXmax','KEXPEP','vemax','KeFBP','ne','acetate','d','kPEPout'};

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
    elseif datatype == 3
        xlabel = parlist(idx(1));
        ylabel = cnclist(idx(2));
    end
    zlabel = {};
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
% plot data from initial value to first bifurcation
% pval = p(ipx,1:xindex(1));
% yval = y(idx,1:xindex(1));
% 
% if any(evalRe(:,xindex(1))<0) && any(evalRe(:,xindex(2))>=0)
%     % entering unstable from stable region 
%     LineP.LineStyle = '-';
%     LineP.Color = 'k';
% elseif any(evalRe(:,xindex(1))>=0) && any(evalRe(:,xindex(2))<0)
%     % entering stable from unstable region
%     LineP.LineStyle = '--';
%     LineP.Color = 'k';
% end
% % plot data
% if ~isempty(hfig)
%     axes(hfig);
%     set(gca,'NextPlot','add');
%     plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
%     hold on
% else
%     plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
%     hold on
% end
% 
% % mssval = zeros(length(xindex),2);
% for ip = 1:length(xindex)
%     if ip<length(xindex)
%         pval = p(ipx,xindex(ip)+1:xindex(ip+1));
%         yval = y(idx,xindex(ip)+1:xindex(ip+1));
%     else
%         pval = p(ipx,xindex(ip):end);
%         yval = y(idx,xindex(ip):end);
%     end
%     
%     % check stability
%      if any(evalRe(:,xindex(ip))<0) && any(evalRe(:,xindex(ip)+1)>=0)
%         % entering unstable from stable region 
%         LineP.LineStyle = '--';
%         LineP.Color = 'k';
%     elseif any(evalRe(:,xindex(ip))>=0) && any(evalRe(:,xindex(ip)+1)<0)
%         % entering stable from unstable region
%         LineP.LineStyle = '-';
%         LineP.Color = 'k';
%      end
%     
% %     if xindex(ip)>1
% %        
% %     else
% %         if any(evalRe(:,xindex(ip))<0) && any(evalRe(:,xindex(ip+1))>=0)
%         % entering unstable from stable region 
%     
% %     mssval(ip,1) = y(idx,xindex(ip));
% %     mssval(ip,2) = p(ipx,xindex(ip));
%     plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
%     line(p(ipx,xindex(ip)),y(idx,xindex(ip)),'LineStyle','none',...
%                              'Marker','o','MarkerEdgeColor','r',...
%                              'MarkerFaceColor','r','MarkerSize',6);
%     hold on
%         
% end



        


