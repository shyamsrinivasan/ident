%Parallel sampling of metabolites and model ensembles
function [ensb,variable] = sampleMet_parallel(FBAmodel,pvector,nmodels)
if nargin < 3
    nmodels = 1;
end
nsmp = 1000;
var = cell(nsmp,1);
Fvar1 = zeros(nsmp,1);
Fvar2 = zeros(nsmp,1);
EScell = cell(nsmp,1);
parfor ismp = 1:nsmp
    [MC,KVl] = sampleMet(FBAmodel);
    [ensb,flag1,flag2] = build_ensemble(1,FBAmodel,pvector,MC,KVl);
%     if ~flag
        %Resample Kms        
%         ensb = resample_ensemble(ensb,FBAmodel,MC);
        var{ismp} = MC;        
        EScell{ismp} = ensb;       
%     end
    Fvar1(ismp) = flag1;
    Fvar2(ismp) = flag2;
end


ensb = struct();
imodel = 1;
for ism = 1:nsmp
    if imodel <= nmodels     
        mname = sprintf('model%d',imodel);
        if Fvar1(ism)==1 && Fvar2(ism)~=1            
            ensb.(mname) = EScell{ism}.model1;
            variable.(mname).MC = var{ism};
            imodel = imodel+1;
        end        
    end
end
