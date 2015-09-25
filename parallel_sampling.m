function mc = parallel_sampling(model,nsample,mc)
if nargin<3
    mc = zeros(model.nt_metab,1);
end
if nargin<2
    nsample = 1;
end

smp = cell(nsample,1);
%extracellular metabolites in M moles/L
met.glc = 0.2;
met.o2 = 1e-5;
met.pi = 1e-6;
met.h = 1e-7;
met.co2 = 1e-8;

for ism = 1:nsample    
    %extracellular metabolites in M moles/L
    smp{ism} = iconcentration(model,met);
    
    %intracellular metabolites for thermodynamic consistency in mM mmoles/L
    smp{ism} = sample_metabolites(model,mc);    
end

if nsample > 1
    mc = struct();
    for ism = 1:nsample
        mid = sprintf('model%s',ism);
        mc.(mid) = smp{ism};
    end
else
    mc = smp{1};
end