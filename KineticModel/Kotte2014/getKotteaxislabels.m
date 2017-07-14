function [xlabel,ylabel,zlabel] = getKotteaxislabels(plotype,datatype,idx)
if nargin <3
    idx = [];
end

fluxlist = {'v1 a.u.','ENZC a.u.','ECbiomass(v3) a.u.',...
           'v2 a.u.','v4 a.u.','ENZR a.u.'};
cnclist = {'pep a.u.','fdp a.u.','ENZ a.u.'};  
parlist = {'K1acetate a.u.','K3fdp a.u.','L3','K3pep a.u.','K2pep a.u.','vemax',...
           'Kefdp a.u.','ne','acetate a.u.','d a.u.','kPEPout a.u.',...
           'k1cat a.u.','v3max a.u.','v2max a.u'};
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
        if idx(2)<=length(cnclist)
            ylabel = cnclist(idx(2));
        else
            ylabel = parlist(idx(2));
        end
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
