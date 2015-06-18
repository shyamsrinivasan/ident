function [vProdLPmax,vProdLPmin,vLPmax,vLPmin,Maxflag,Minflag] =...
         solveLP(model,prod,enz,bounds,prxnid)
%      if nargin < 5
%          if ~isempty(prod)
%              prdrxn = model.S(strcmpi(prod,model.Metabolites),:)
%                 
%              prxnid = strcmpi(
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
nint_metab = model.nint_metab-length(find(model.S(:,model.bmrxn)>0));
S = model.S;%(1:nint_metab,:);
%Add additional reactions for reversible reactions
% if any(model.reversible)
%     revCol = S(:,logical(model.reversible));
%     S = [S -revCol];
% end

%Add elements to S for growth dilution
if prxnid ~= model.bmrxn
    bm_col = S(:,model.bmrxn);
    bm_col(bm_col < 0) = bm_col(bm_col < 0 )-1;
    bm_col(bm_col == 0) = -1;
    S(:,model.bmrxn) = bm_col;
end

[nm,nr] = size(S);

%Uptake Fluxes
if isfield(bounds,'Vuptake')
    if ~isempty(bounds.Vuptake)
        vl(model.Vupind) = bounds.Vuptake;
        vu(model.Vupind) = bounds.Vuptake;
    end
end

Mglc = strcmpi('glc[e]',model.Metabolites);
Vglc = find(model.S(Mglc,:)<0);
vl(Vglc) = bounds.Vuptake;
vu(Vglc) = bounds.Vuptake;

%Growth Fluxes
if prxnid ~= model.bmrxn
    if isfield(model,'gmax')
        vl(model.bmrxn) = model.gmax;
        vu(model.bmrxn) = model.gmax;
    end
else
    vl(model.bmrxn) = 0;
end
%Exchnage Reactions Cannot occur in reverse
% vl(model.Vex) = 0;

%Default Bounds - Reversible reactions

b = zeros(nm,1);
%Exchange of Product
% vl(model.Vexind) = 0;
%Objective Function
cprod = sparse(1,prxnid,1,1,nr);

[vLPmax,vProdLPmax,Maxflag] = cplexlp(-cprod(:),[],[],S,b,vl,vu);
[vLPmin,vProdLPmin,Minflag] = cplexlp(cprod(:),[],[],S,b,vl,vu);

return

