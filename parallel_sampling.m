function [mc,pvec,smp] = parallel_sampling(model,pvec,nsample)
if nargin<3
    nsample = 500;
end

%generate one metabolite concentration for parameter estimation
%get one set of concentrations and coresponding delGr
[mc,assignFlag,delGr,vCorrectFlag] = getiConEstimate(model);

%extracellular metabolites in M moles/L
met.glc = 0.2;
met.o2 = 1e-5;
met.pi = 1e-3;
met.h = 1e-7;
met.co2 = 1e-8;
met.h2o = 55.0;

mc = iconcentration(model,met,mc,assignFlag);
pvec.delGr = delGr;

%sample met using ACHR
%get more than one set of concentrations using ACHR sampling
clear assignFlag
pts = [];
ptsdelGr = [];
% [pts,assignFlag,ptsdelGr,vCorrectFlag] =...
% ACHRmetSampling(model,1,nsample,200);

% pts = iconcentration(model,met,pts,assignFlag);

if ~isempty(pts)
    smp = cell(nsample,1);
    for ism = 1:nsample    
        %extracellular metabolites in M moles/L

        smp{ism,1} = pts(:,ism);
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
