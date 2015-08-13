function [vProdLPmax,vProdLPmin,vLPmax,vLPmin,Maxflag,Minflag] =...
         solveLP(model,prod,enz,bounds,prxnid,fixgrowth)
if nargin < 6
    fixgrowth = 0;
end
if nargin < 4
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

%Uptake Fluxes
if isfield(bounds,'Vuptake')
    if ~isempty(bounds.Vuptake)
        vl(model.VFup) = bounds.Vuptake;
        vu(model.VFup) = bounds.Vuptake;
    end
end
%Set uptake to ireversible
if ~isempty(model.Vupind)
    vl(model.Vupind) = 0;
    vu(model.Vupind) = 100;
end
%Set VFex to irreversible
if ~isempty(model.VFex)
    vl(model.VFex) = 0;
    vu(model.VFex) = 100;
end
%Set Vex to irreversible
if ~isempty(model.VFex)
    vl(model.Vex) = 0;
    vu(model.Vex) = 100;
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
%Exchnage Reactions Cannot occur in reverse
% rxnid = strcmpi('PPS',model.rxns);
% vl(rxnid) = 0;
% vu(rxnid) = 0;
%Default Bounds - Reversible reactions

b = zeros(nm,1);
%Exchange of Product
% vl(model.Vexind) = 0;
%Objective Function
% prxnid = find(strcmpi('EXpep',model.rxns));
% pps = strcmpi('PPS',model.rxns);
% vl(pps) = 0;
% vu(pps) = 0;

cprod = sparse(1,prxnid,1,1,nr);


[vLPmax,vProdLPmax,Maxflag] = cplexlp(-cprod(:),[],[],S,b,vl,vu);
if Maxflag > 0
    if abs(vProdLPmax)<1e-4
        vProdLPmax = 0;
    end
end
[vLPmin,vProdLPmin,Minflag] = cplexlp(cprod(:),[],[],S,b,vl,vu);

return

