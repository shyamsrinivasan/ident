%make multiple perturbations in parallel
function multiplePerturbations(idx,val,initval)
if nargin<3
    error('mulPert:NoInitVal','No initial value given for perturbation');
end

np = length(idx);
if np~=length(val)
    error('mulPert:NoVal','Length of value and index vectors do not match');
end

for id = 1:np
    Nimc = perturbEqSolution(model,initval,[],[],idx(id),val(id));
    