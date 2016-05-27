function FIGmssdynamicswpvec(y,tout,pvec,idx,idp,jpts,type)
% function to plot dynamic time profiles of flux/concentration for systems
% whose multistationarity has been evaluated using MATCONT
switch type
    case 'flux'
        type = 1;
    case 'conc'
        type = 2;
    otherwise
        type = 0;
end
if length(size(y))>2
    % dynamic data
else
    % ss data
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
if ndp > 1
    ncol = 2;
else
    ncol=1;
end

% sort data as per increasing pvec

ifig = 1;
for ip = 1:ndp
    [~,index] = sortrows(pvec,idp(ip));
    id = ismember(index,jpts);
    select_sorted_p = pvec(index(id),:);
    pval = select_sorted_p(:,idp(ip));
    for iy = 1:ndx
        ydata = reshape(y(idx(iy),:,index(id),ip),length(tout),length(jpts));
        pvec_select = pvec(index(id),:);
        xdata = repmat(tout,1,length(jpts));
        hsfig = subplot(nrows,ncol,ifig);
        set(hsfig,'NextPlot','add');
        ptname = cell(length(jpts),1);
        for ipts = 1:length(jpts)   
            ptname{ipts} = [num2str(jpts(ipts)) '=' num2str(pval(ipts))];             
        end           
        ht = line(xdata,ydata);
        text(xdata(end,:),ydata(end,:),ptname);
        set(ht,{'DisplayName'},ptname);
        
        if type == 1
            switch idx(iy)
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
        elseif type == 2
            switch idx(iy)
                case 1
                    ylabel = sprintf('PEP');
                case 2
                    ylabel = sprintf('FDP');
                case 3
                    ylabel = sprintf('ENZ');                
            end 
        end
        xlabel = sprintf('Time (s)');
        set(get(gca,'YLabel'),'String',ylabel);
        set(get(gca,'XLabel'),'String',xlabel);
        ifig = ifig+1;
    end    
end
