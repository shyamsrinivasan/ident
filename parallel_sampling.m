function [mc,pvec,smp] = parallel_sampling(model,pvec,nsample)
if nargin<3
    nsample = 500;
end

%extracellular metabolites in M moles/L
% met.glc = 0.2;
met.pi_e = 1e-2;
met.o2_e = 0.0025;
met.h_e = 1e-7;
met.co2_e = 0.002;%1e-8;
met.h2o_e = 55.0;%55.0;
met.pi_c = 1e-3;
met.h_c = 1e-7;
met.h2o_c = 55.0;%1.53e-13;

fprintf('Generating single feasible concentration sample\n');
%generate one metabolite concentration for parameter estimation
%get one set of concentrations and coresponding delGr
[mc,assignFlag,delGr,model,vCorrectFlag] = getiConEstimate(model);

[mc,assignFlag] = iconcentration(model,met,mc,assignFlag);
%M to mM
mc = mc*1000;
pvec.delGr = delGr;

% fprintf('Generating %d concentration samples using ACHR\n',nsample);
%sample met using ACHR
%get more than one set of concentrations using ACHR sampling
clear assignFlag
pts = [];
ptsdelGr = [];
% [pts,assignFlag,ptsdelGr,vCorrectFlag] =...
% ACHRmetSampling(model,1,nsample,200);

% pts = iconcentration(model,met,pts,assignFlag);

model.lb(strcmpi(model.mets,'pi[e]')) = 1e-9;
model.lb(strcmpi(model.mets,'o2[e]')) = 1e-9;
model.lb(strcmpi(model.mets,'co2[e]')) = 1e-9;
model.lb(strcmpi(model.mets,'h[e]')) = 1e-9;
model.lb(strcmpi(model.mets,'h2o[e]')) = 1e-9;

model.ub(strcmpi(model.mets,'pi[e]')) = 250;
model.ub(strcmpi(model.mets,'o2[e]')) = 250;
model.ub(strcmpi(model.mets,'co2[e]')) = 250;
model.ub(strcmpi(model.mets,'h[e]')) = 250;
model.ub(strcmpi(model.mets,'h2o[e]')) = 250;

if ~isempty(pts)
    smp = cell(nsample,1);
    for ism = 1:nsample    
        %extracellular metabolites in M moles/L
        %from M to mM
        smp{ism,1} = pts(:,ism)*1000;
        smp{ism,2} = pvec;
        smp{ism,2}.delGr = ptsdelGr(:,ism);

        %details in structure format
    %     if nsample>1
    %         mid = sprintf('model%s',ism);
    %         variable.(mid) = smp{ism,1};
    %         pvec.(mid).delGr = ptsdelGr(:,ism); 
    %     else
    %         variable = mc;
    %         pvec.delGr = delGr;
    %     end
    end
else
    smp = {};
end

% if nsample > 1
%     mc = struct();
%     newpvec = struct();    
%     for ism = 1:nsample
%         mid = sprintf('model%s',ism);
%         mc.(mid) = smp{ism};
%         newpvec.(mid) = pvec;
%         newpvec.(mid).delGr = smp{ism,2};
%     end
% else
%     mc = smp{1,1};
%     pvec.delGr = smp{1,2};
% end
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\ecoliN1_MC1.mat');
% load('C:\Users\shyam\Documents\Courses\CHE1125Project\mat_files\KineticModel\ecoliN1_pvec1.mat');
