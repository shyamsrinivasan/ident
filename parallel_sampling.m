function [mc,pvec,smp] = parallel_sampling(model,pvec,funName,met,mc,rxn_add,nsample)
if nargin<7
    nsample = 500;
end
if nargin<6
    rxn_add = {};
end
if nargin<5
    mc = [];
end
if nargin<4 || isempty(met)
    met = struct([]);
end


if nargin == 6
    % generate one metabolite concentration for parameter estimation
    % get one set of concentrations and coresponding delGr
    [mc,assignFlag,delGr,model,vCorrectFlag] = getiConEstimate(model,funName,mc,rxn_add);
    [mc,assignFlag] = iconcentration(model,met,mc,assignFlag);
    %M to mM
    mc = mc*1000;
    pvec.delGr = delGr;
    if nargout == 3
        smp = {};
    end
elseif nargin>6
    % sample met using ACHR
    % get more than one set of concentrations using ACHR sampling
    [pts,assignFlag,ptsdelGr] =...
    ACHRmetSampling(model,funName,mc,rxn_add,1,nsample,200);
    pts = iconcentration(model,met,pts,assignFlag);
    if ~isempty(pts)
        smp = cell(size(pts,2),1);
        for ism = 1:size(pts,2)
            %extracellular metabolites in M moles/L
            %from M to mM
            smp{ism,1} = pts(:,ism)*1000;
            smp{ism,2} = pvec;
            smp{ism,2}.delGr = ptsdelGr(:,ism);        
        end        
    else
        smp = {};
    end
    if nargout==1
        mc = smp;
    elseif nargout==3
        mc = [];
        pvec = [];
    end
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
