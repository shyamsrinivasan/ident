function new_mc = perturbEqSolution(model,mc,change_pos,change_neg,idx,val)
if nargin<6
    val = 10;
end
if nargin<5
    idx = [];
end
%decrease in concentration
if nargin<4
    change_neg = [];
end
%increase in concentration
if nargin<3
    change_pos = [];
end

if ~isempty(idx)
    for id = 1:length(idx)
        if val(id)>0
            change_pos.(model.mets{idx(id)}) = val(id);
        elseif val(id)<0
            change_neg.(model.mets{idx(id)}) = val(id);
        end
    end
end
%backup mc
old_mc = mc;

%concentration increase
if ~isempty(change_pos)
     new_mc = changeInitialCondition(model,old_mc,change_pos);
end
%concentration decrease
if ~isempty(change_neg)
    new_mc = changeInitialCondition(mdoel,old_mc,[],change_neg);
end

if isempty(change_pos) && isempty(change_neg) && isempty(idx) 
    allMets = 1;
    new_mc = changeInitialCondition(model,old_mc,[],[],allMets);
end
    
   


 