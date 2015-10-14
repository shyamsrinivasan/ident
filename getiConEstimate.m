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
% if size(bounds.A,2) == length(bounds.mets)
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
    
    %solve once more to  obtain consistent concentrations  size(bounds.A,2)  
    LPmax = solvemetLP(bounds);
    if LPmax.flag>0 
        %do not include slack variables
        [mc,assignFlag] = separate_slack(LPmax.x,model,bounds);
        mc = exp(mc);
    else
        error('mcEst:LPinfeas',...
            'LP for thermodynamic metabolite conentrations is infeasible');
    end
%     if ~isempty(x)
%         [mc,assignFlag] = assignConc(x,model,bounds);        
%     end
    if ~isempty(xmax)
%         mcMax = assignConc(xmax,model,bounds);
        mcMax = separate_slack(xmax,model,bounds);
        mcMax = exp(mcMax);
    end
    if ~isempty(xmin)
%         mcMin = assignConc(xmin,model,bounds);
        mcMin = separate_slack(xmin,model,bounds);
        mcMin = exp(mcMin);
    end
    %check for delGr values
    delGr = checkdelGr(model,mc);    
    RxnDir = delGr.*model.Vss;
% else
%     error('getiConEst:sizeCheck',...
%         'Number of metabolites in bounds.S and bounds.mets do not match');
% end

%concentrations for thermodynamically inactive reactions
%delG = 0 and Vss = 0


% function [mc,assignFlag] = assignConc(x,model,bounds)
% %assign concentrations to the model at large
% mc = zeros(length(model.mets),1);
% assignFlag = zeros(length(model.mets),1);
% for im = 1:length(bounds.mets)
%     tfm = strcmpi(bounds.mets{im},model.mets);
%     if any(tfm)
%         mc(tfm) = x(im);
%         assignFlag(tfm) = 1;
%     end
% end

function delGr = checkdelGr(model,mc)

%eliminate consideration for excess cofators
%pi[c],pi[e],h[c],h[e],h2o[c]
pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));

vspl = [find(strcmpi(model.rxns,'THD2'))...
        find(strcmpi(model.rxns,'NADH16'))...
        find(strcmpi(model.rxns,'ATPS4r'))];

%check for delGr values
Vind = model.Vind;
delGr = zeros(model.nt_rxn,1);
for irxn = 1:length(model.Vind)
    if ismember(Vind(irxn),vspl)
        q8 = find(strcmpi(model.mets,'q8[c]'));
        q8h2 = find(strcmpi(model.mets,'q8h2[c]'));
    else
        q8 = [];
        q8h2 = [];
    end
    sbid = model.S(:,Vind(irxn))<0;
    prid = model.S(:,Vind(irxn))>0;
    sbid([pic pie hc he h2o q8 q8h2]) = 0;
    prid([pic pie hc he h2o q8 q8h2]) = 0;
    sb = prod(mc(sbid));
    pr = prod(mc(prid));
    delGr(Vind(irxn)) = 0.008314*298*log(pr/(sb*model.Keq(Vind(irxn))));
    fprintf('delG = %2.3g \t Vss = %3.4g \n',delGr(Vind(irxn)),model.Vss(Vind(irxn)));
    
end
