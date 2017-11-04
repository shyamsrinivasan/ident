function [model,bounds] = changebounds_toy(model,ess_rxn,bounds,fixgrowth)
if nargin<4
    fixgrowth = 0;
end
if nargin<3
    bounds = struct();
    if isfield(model,'Vuptake')
        Vuptake = model.Vuptake;
    end
    nr = size(model.S,2);
    vl = zeros(nr,1);
    vl(logical(model.rev)) = -100;
    vu = zeros(nr,1);
    vu(vu==0);
else
    if isfield(bounds,'Vuptake')
        Vuptake = bounds.Vuptake;
    end
    vl = bounds.vl;
    vu = bounds.vu;
end
if nargin<2
    ess_rxn = {};
end

%Uptake Fluxes
if any(Vuptake)
    vl(logical(Vuptake)) = -Vuptake(logical(Vuptake));
end

%change bounds for exchange metabolites
essid = [];
for iess = 1:length(ess_rxn)
    essid = union(essid,find(strcmpi(ess_rxn{iess},model.rxns)));
end
if isfield(model,'VFex')
    Vess = setdiff(model.VFex,essid);
    vl(Vess) = 0;
end

bounds.Vuptake = Vuptake;
bounds.vl = vl;
bounds.vu = vu;

