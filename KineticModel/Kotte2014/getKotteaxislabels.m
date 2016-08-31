function [xlabel,ylabel,zlabel] = getKotteaxislabels(plotype,datatype,idx)
if nargin <3
    idx = [];
end

fluxlist = {'ACpts a.u','ENZC a.u','ECbiomass(FDP) a.u',...
           'GLUX a.u','PEPout a.u'};
cnclist = {'PEP a.u','FDP a.u','ENZ a.u'};  
parlist = {'KEacetate','KFbpFBP','Lfbp','KFbpPEP','KEXPEP','vemax',...
           'KeFBP','ne','acetate a.u','d a.u','kPEPout a.u',...
           'kEcat a.u','vFbpmax a.u','vEXmax a.u'};
% parlist = {'Parameter 1, kEcat','','','Parameter 2, vFbpmax','','',...
%            'Parameter 3, vEXmax','','','','','','','Parameter 4, kPEPout'};

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
elseif plotype == 2
    if datatype == 1
        xlabel = fluxlist(idx(1));
        ylabel = fluxlist(idx(2));
%         zlabel = fluxlist(idx);
%         ylabel = parlist(idp(2));
    elseif datatype == 2
        xlabel = cnclist(idx(1));
        ylabel = cnclist(idx(2));
    elseif datatype == 3
%         xlabel = parlist(idx(1)-length(cnclist));
        xlabel = parlist(idx(3));
        ylabel = cnclist(idx(2));
    elseif datatype == 4
%         xlabel = parlist(idx(1)-length(fluxlist));
        xlabel = parlist(idx(3));
        ylabel = fluxlist(idx(2));    
    end
    zlabel = {};
end
