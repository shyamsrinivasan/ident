function [xlabel,ylabel,zlabel] = getaxislabels(plotype,datatype,model,idx)
if nargin<4
    idx = [];
end

flxlst = model.rxns;
metlst = model.mets;

if plotype==2
    switch datatype
        case 11
            xlabel = {'Time s'};
            ylabel = metlst(idx(2));
        case 12
            xlabel = {'Time s'};
            ylabel = flxlst(idx(2));
        case 1
            xlabel = metlst(idx(1));
            ylabel = metlst(idx(2));
        case 2
            xlabel = flxlst(idx(1));
            ylabel = flxlst(idx(2));
        case 3
%             xlabel = parlist(idx(3));
%             ylabel = metlst(idx(2));
        case 4
%             xlabel = parlist(idx(3));
%             ylabel = flxlst(idx(2));   
    end    
    zlabel = {};
elseif plotype==3
end
