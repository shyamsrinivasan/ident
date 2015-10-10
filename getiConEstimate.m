function [mc,assignFlag,mcMax,mcMin] = getiConEstimate(model)

%setup problem with default constraints for thermodynamically active
%reaction
% delG > or < 0 and Vss < or > 0
bounds = setupMetLP(model);
x=[];
xmax = [];
xmin = [];
mcMax = [];
mcMin = [];

%check
if ~isfield(bounds,'A')
    error('getiest:NoA','No stoichiometric transpose found');
end
if size(bounds.A,2) == length(bounds.mets)
%     nmets = size(bounds.A,1);    
    
%     %loop through all metabolites to get min and max
%     xmax = zeros(nmets,1);
%     xmin = zeros(nmets,1);
%     for im = 1:nmets
%         prxnid = im;
%         [LPmax,LPmin] = solvemetLP(bounds,prxnid);
%         if LPmax.flag>0
%             xmax(im) = exp(LPmax.obj);
%         end
%         if LPmin.flag>0
%             xmin(im) = exp(LPmin.obj);
%         end
%     end
    
    %solve once more to  obtain consistent concentrations    
    LPmax = solvemetLP(bounds);
    if LPmax.flag>0 
        x = exp(LPmax.x);
    else
        error('mcEst:LPinfeas',...
            'LP for thermodynamic metabolite conentrations is infeasible');
    end
    if ~isempty(x)
        [mc,assignFlag] = assignConc(x,model,bounds);        
    end
    if ~isempty(xmax)
        mcMax = assignConc(xmax,model,bounds);
    end
    if ~isempty(xmin)
        mcMin = assignConc(xmin,model,bounds);
    end
else
    error('getiConEst:sizeCheck',...
        'Number of metabolites in bounds.S and bounds.mets do not match');
end

%concentrations for thermodynamically inactive reactions
%delG = 0 and Vss = 0


function [mc,assignFlag] = assignConc(x,model,bounds)
%assign concentrations to the model at large
mc = zeros(length(model.mets),1);
assignFlag = zeros(length(model.mets),1);
for im = 1:length(bounds.mets)
    tfm = strcmpi(bounds.mets{im},model.mets);
    if any(tfm)
        mc(tfm) = x(im);
        assignFlag(tfm) = 1;
    end
end