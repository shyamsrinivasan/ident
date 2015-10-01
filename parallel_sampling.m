function [mc,pvec] = parallel_sampling(model,pvec,nsample)
if nargin<3
    nsample = 1;
end

smp = cell(nsample,1);
%extracellular metabolites in M moles/L
met.glc = 0.2;
met.o2 = 1e-5;
met.pi = 1e-3;
met.h = 1e-7;
met.co2 = 1e-8;
met.h2o = 55.0;

for ism = 1:nsample    
    %extracellular metabolites in M moles/L
    smp{ism} = iconcentration(model,met);
    
    %intracellular metabolites for thermodynamic consistency in mM mmoles/L
    [mc,delGr] = metsample(model,smp{ism});
    smp{ism,1} = mc;
    smp{ism,2} = delGr;
end

if nsample > 1
    mc = struct();
    newpvec = struct();    
    for ism = 1:nsample
        mid = sprintf('model%s',ism);
        mc.(mid) = smp{ism};
        newpvec.(mid) = pvec;
        newpvec.(mid).delGr = smp{ism,2};
    end
else
    mc = smp{1,1};
    pvec.delGr = smp{1,2};
end
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\ecoliN1_MC1.mat');
load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\ecoliN1_pvec1.mat');
