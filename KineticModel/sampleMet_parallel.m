%Parallel sampling of metabolites and model ensembles
function [var,Fvar,EScell] = sampleMet_parallel(FBAmodel,parameter,variable,nmodels)
if nargin < 4
    nmodels = 1;
end
nsmp = 500000;
var = cell(nsmp,1);
Fvar = zeros(nsmp,1);
EScell = cell(nsmp,1);
for ismp = 1:nsmp
    MC = sampleMet(FBAmodel);
    [ensb,flag] = build_ensemble(nmodels,FBAmodel,parameter,MC);
%     if ensb.model1.Vmax(9) > 0 && ensb.model1.Vmax(7)>0
%         fprintf('Yay! %d',ismp);
%     end
    if ~flag
        %Resample Kms        
        ensb = resample_ensemble(ensb,FBAmodel,MC);
        var{ismp} = MC;        
        EScell{ismp} = ensb;       
    end
    Fvar(ismp) = flag;
end

for ism = 1:nsmp
    if Fvar(ism)~=1
        ensb = EScell{ism};
        variable.MC = var{ism};
    end
end
