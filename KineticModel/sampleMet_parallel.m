%Parallel sampling of metabolites and model ensembles
function [ensb,variable] = sampleMet_parallel(FBAmodel,parameter,nmodels)
if nargin < 3
    nmodels = 1;
end
nsmp = 100;
var = cell(nsmp,1);
Fvar = zeros(nsmp,1);
EScell = cell(nsmp,1);
parfor ismp = 1:nsmp
    MC = sampleMet(FBAmodel);
    [ensb,flag] = build_ensemble(1,FBAmodel,parameter,MC);
    if ~flag
        %Resample Kms        
        ensb = resample_ensemble(ensb,FBAmodel,MC);
        var{ismp} = MC;        
        EScell{ismp} = ensb;       
    end
    Fvar(ismp) = flag;
end


ensb = struct();
imodel = 1;
for ism = 1:nsmp
    if imodel <= nmodels     
        mname = sprintf('model%d',imodel);
        if Fvar(ism)~=1            
            ensb.(mname) = EScell{ism}.model1;
            variable.(mname).MC = var{ism};
            imodel = imodel+1;
        end        
    end
end
