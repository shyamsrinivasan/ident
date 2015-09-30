function [vProdLPmax,vLPmax,vProdLPmin,vLPmin,Maxflag,Minflag,model] =...
         solveLP(model,bounds,prxnid,fixgrowth)
if nargin < 4
    fixgrowth = 0;
end
if nargin < 2
    nr = size(model.S,2);
    vl = zeros(nr,1);
    vl(vl==0) = -100;
    vu = zeros(nr,1);
    vu(vu==0) = 100;
else
    vl = bounds.vl;
    vu = bounds.vu;
end
% nint_metab = model.nint_metab-length(find(model.S(:,model.bmrxn)>0));
S = model.S;%(1:nint_metab,:);
%Add additional reactions for reversible reactions
% if any(model.reversible)
%     revCol = S(:,logical(model.reversible));
%     S = [S -revCol];
% end

%Add elements to S for growth dilution
% bm_col = S(:,model.bmrxn);
% bm_col(bm_col < 0) = bm_col(bm_col < 0 )-1;
% bm_col(bm_col == 0) = -1;
% S(:,model.bmrxn) = bm_col;

[nm,nr] = size(S);
%Flux calculation based on intial concentrations
% flux = ExFlux(model,Y,zeros(nr,1),model.Vupind,'mm');

%change bounds for exchange metabolites
ess_rxn = {'exCO2','exH','exH2O','exPI','exO2'};
essid = [];
for iess = 1:length(ess_rxn)
    essid = union(essid,find(strcmpi(ess_rxn{iess},model.rxns)));
end
Vess = setdiff(model.VFex,essid);
vl(Vess) = 0;

%atp maintanance
vl(strcmpi(model.rxns,'ATPM')) = 8.39;

%Uptake Fluxes
if isfield(bounds,'Vuptake')
    if ~isempty(bounds.Vuptake)
        vl(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
        vl(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
        vu(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
        vu(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
    end
end

%Growth Fluxes
if fixgrowth
    if prxnid ~= model.bmrxn
        if isfield(model,'gmax')
            vl(model.bmrxn) = model.gmax;
            vu(model.bmrxn) = model.gmax;
        end
    end
else
    vl(model.bmrxn) = 0;    
end

b = zeros(nm,1);
%Exchange of Produ ct
% vl(model.Vexind) = 0;
%Objective Function
% prxnid = find(strcmpi('EXpep',model.rxns));
% pps = strcmpi('PPS',model.rxns);
% vl(pps) = 0;
% vu(pps) = 0;

cprod = sparse(1,prxnid,1,1,nr);
% cprod = sparse(1,73,1,1,nr);

[vLPmax,vProdLPmax,Maxflag] = cplexlp(-cprod(:),[],[],S,b,vl,vu);
if Maxflag > 0
    if abs(vProdLPmax)<1e-4
        vProdLPmax = 0;
    end
end
[vLPmin,vProdLPmin,Minflag] = cplexlp(cprod(:),[],[],S,b,vl,vu);
model.lb = vl;
model.ub = vu;
model.c = cprod;
model.b = b;

return

