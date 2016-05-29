function plotPVECvsX(pvec,idp,x,idx)

figure
% nrxn = size(x,1);
ndp = length(idp);
ndx = length(idx);
nfig = ndx*ndp;
if rem(nfig,2)~=0
    nfig = nfig+1;
end

nrows = nfig/2;
if ndp > 1
    ncol = 2;
else
    ncol=1;
end

ifig = 1;
for ip = 1:ndp
    for ix = 1:ndx
        hsubfig = subplot(nrows,ncol,ifig);
        line(pvec(:,idp(ip)),x(idx(ix),:),'LineStyle','none',...
                         'Marker','o','MarkerEdgeColor','k',...
                         'MarkerFaceColor','k','MarkerSize',7);
        switch idx(ix)
            case 1
                ylabel = sprintf('ACpts mmole/h');
            case 2
                ylabel = sprintf('ENZC mmole/h');
            case 3
                ylabel = sprintf('ECbiomass(FDP) mmole/h');
            case 4
                ylabel = sprintf('GLUX mmole/h');
            case 5
                ylabel = sprintf('PEPout mmole/h');
        end    
        ptid = ['P ' num2str(ip)];
        xlabel = sprintf('Parameter Values %s',ptid);
        set(get(gca,'YLabel'),'String',ylabel);
        set(get(gca,'XLabel'),'String',xlabel);
        ifig = ifig + 1;
    end
end




% while ifig <= ndp
%     hsubfig = subplot(nrows,ncol,ifig);    
%     line(pvec(),x(idp(ifig),:),'LineStyle','none',...
%                              'Marker','o','MarkerEdgeColor','k',...
%                              'MarkerFaceColor','k','MarkerSize',7);
%     switch idp(ifig)
%         case 1
%             ylabel = sprintf('ACpts mmole/h');
%         case 2
%             ylabel = sprintf('ENZC mmole/h');
%         case 3
%             ylabel = sprintf('ECbiomass(FDP) mmole/h');
%         case 4
%             ylabel = sprintf('GLUX mmole/h');
%         case 5
%             ylabel = sprintf('PEPout mmole/h');
%     end        
%     xlabel = sprintf('Parameter Values');
%     set(get(gca,'YLabel'),'String',ylabel);
%     set(get(gca,'XLabel'),'String',xlabel);
% %     ht = text(
%     ifig = ifig+1;    
% end