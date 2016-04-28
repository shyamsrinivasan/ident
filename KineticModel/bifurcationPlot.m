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
    xindex = xindex(2:end-1);
end
evalRe = real(f1);

% plot data from initial value to first bifurcation
pval = p(ipx,1:xindex(1));
yval = y(idx,1:xindex(1));
if any(evalRe(:,xindex(1))<0) && any(evalRe(:,xindex(2))>=0)
    % entering unstable from stable region 
    LineP.LineStyle = '-';
    LineP.Color = 'k';
elseif any(evalRe(:,xindex(1))>=0) && any(evalRe(:,xindex(2))<0)
    % entering stable from unstable region
    LineP.LineStyle = '--';
    LineP.Color = 'k';
end
% plot data
if ~isempty(hfig)
    axes(hfig);
    set(gca,'NextPlot','add');
    plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
    hold on
else
    plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
    hold on
end

for ip = 1:length(xindex)
    if ip<length(xindex)
        pval = p(ipx,xindex(ip)+1:xindex(ip+1));
        yval = y(idx,xindex(ip)+1:xindex(ip+1));
    else
        pval = p(ipx,xindex(ip):end);
        yval = y(idx,xindex(ip):end);
    end
    
    % check stability
     if any(evalRe(:,xindex(ip))<0) && any(evalRe(:,xindex(ip)+1)>=0)
        % entering unstable from stable region 
        LineP.LineStyle = '--';
        LineP.Color = 'k';
    elseif any(evalRe(:,xindex(ip))>=0) && any(evalRe(:,xindex(ip)+1)<0)
        % entering stable from unstable region
        LineP.LineStyle = '-';
        LineP.Color = 'k';
     end
    
%     if xindex(ip)>1
%        
%     else
%         if any(evalRe(:,xindex(ip))<0) && any(evalRe(:,xindex(ip+1))>=0)
        % entering unstable from stable region 
    
    plot(pval,yval,'LineStyle',LineP.LineStyle,'Color',LineP.Color,'LineWidth',3);
    line(p(ipx,xindex(ip)),y(idx,xindex(ip)),'LineStyle','none',...
                             'Marker','o','MarkerEdgeColor','r',...
                             'MarkerFaceColor','r','MarkerSize',6);
    hold on
        
end



        


