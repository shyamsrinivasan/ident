function bifurcationPlot(y,p,s1,f1,idx,ipx,hfig)
if nargin<7
    hfig = [];
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
    if (kndex+1)<=length(xindex)
        pval = p(ipx,xindex(kndex):xindex(kndex+1));
        yval = y(idx,xindex(kndex):xindex(kndex+1));        
    else
        pval = p(ipx,xindex(kndex):end);
        yval = y(idx,xindex(kndex):end);
    end
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
    
    % plot data
    if ~isempty(hfig)
        figure(hfig);
        set(gca,'NextPlot','add');
        plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
        hold on
    else
        plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
        hold on
    end
end

% plot bifurcation points
xindex = xindex(2:end-1);
for ip = 1:length(xindex)
    line(p(ipx,xindex(ip)),y(idx,xindex(ip)),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
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



        


