function [vLPmax,vLPmin,model] =...
         solveLP(model,bounds,prxnid,fixgrowth)
if nargin < 4
    fixgrowth = 0;
end
if nargin < 2
    nr = size(model.S,2);
    if ~isfield(model,'vl')
        vl = zeros(nr,1);
        vl(vl==0) = -100;
    else
        vl = model.vl;
    end
    if ~isfield(model,'vu')
        vu = zeros(nr,1);
        vu(vu==0) = 100;
    else
        vu = model.vu;
    end
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

[nm,nr] = size(S);
%Flux calculation based on intial concentrations
% flux = ExFlux(model,Y,zeros(nr,1),model.Vupind,'mm');

%change bounds for exchange metabolites
ess_rxn = {'exCO2','exH','exH2O','exPI','exO2','exGLC'};
essid = [];
for iess = 1:length(ess_rxn)
    essid = union(essid,find(strcmpi(ess_rxn{iess},model.rxns)));
end
if isfield(model,'VFex')
    Vess = setdiff(model.VFex,essid);
    vl(Vess) = 0;
end

%atp maintanance
vl(strcmpi(model.rxns,'ATPM')) = 8.39;
vu(strcmpi(model.rxns,'ATPM')) = 100;

%ETC
% vl(strcmpi(model.rxns,'ATPS4r')) = -100;
% vu(strcmpi(model.rxns,'ATPS4r')) = 0;
vl(strcmpi(model.rxns,'NADTRHD')) = 0;
vu(strcmpi(model.rxns,'NADTRHD')) = 100;
model.rev(strcmpi(model.rxns,'NADTRHD')) = 0;
vl(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'SUCCt2_2')) = 0;
vu(strcmpi(model.rxns,'FORt2')) = 0;

%Uptake Fluxes
if isfield(bounds,'Vuptake')
    if ~isempty(bounds.Vuptake)
        vl(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
        vl(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
%         vu(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
%         vu(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
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

if ~isfield(model,'b')
    b = zeros(nm,1);
else
    b = model.b;
end

if ~isfield(model,'vl')
    model.vl = vl;
end
if ~isfield(model,'vu')
    model.vu = vu;
end
%Exchange of Produ ct
% vl(model.Vexind) = 0;
%Objective Function
% prxnid = find(strcmpi('EXpep',model.rxns));
% pps = strcmpi('PPS',model.rxns);
% vl(pps) = 0;
% vu(pps) = 0;

cprod = sparse(1,prxnid,1,1,nr);

%maximization
[vmax,vobj,Maxflag] = cplexlp(-cprod(:),[],[],S,b,vl,vu);
if Maxflag > 0
    if abs(vobj)<1e-4
        vobj = 0;
    end
    vLPmax.obj = vobj;
    vLPmax.v = vmax;    
else
    vLPmax.obj = [];
    vLPmax.v = [];    
end
vLPmax.flag = Maxflag;

%minimization
[vmin,vobj,Minflag] = cplexlp(cprod(:),[],[],S,b,vl,vu);
if Minflag > 0
    if abs(vobj)<1e-4
        vobj = 0;
    end
    vLPmin.obj = vobj;
    vLPmin.v = vmin;
else
    vLPmin.obj = [];
    vLPmin.v = [];    
end
vLPmin.flag = Maxflag;

model.vl = vl;
model.vu = vu;
model.c = cprod;
model.b = b;

return

