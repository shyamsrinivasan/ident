function [newmodel,newpvec,newmc,cnstmet] = remove_eMets(model,pvec,mc,rxnid,metid)
% Inputs
% rxnid     rxns that need to be kept in the new model (do not remove these
%           reactions whatsoever
% metid     metabolites whose concentrations are to fixed and removed as
%           variables from the model

if nargout>=4
    cnstmet = [];
end

newmodel = model;
newpvec = pvec;
newmc = mc;

% only keep rxns in rxnid
if ~isempty(rxnid)
    rmvlst = model.rxns(setdiff(1:length(model.rxns),rxnid));
    [newmodel,newpvec,newmc] = addremoveRxns(newmodel,rmvlst,newpvec,newmc);
end


% set mets in metid as constant parameters
if ~isempty(metid)
    cnstmet = zeros(length(newmodel.mets),1);    
    for im = 1:length(metid)
        tfm = strcmpi(newmodel.mets,metid{im});
        if any(tfm)
            cnstmet(tfm) = 1;            
        end
    end    
end

% rearrange metabolites such that variables are first and constants come
% at the end
cnstmet = logical(cnstmet);
varmet = ~cnstmet;
bm_ind = [];
newmodel.mets = [newmodel.mets(varmet,1);...
                 newmodel.mets(cnstmet,1)];   
newmodel.S = [newmodel.S(varmet,:);newmodel.S(cnstmet,:)];
newmodel.SI = [newmodel.SI(varmet,:);newmodel.SI(cnstmet,:)];

if isfield(newmodel,'CMPS')
    newmodel.CMPS = [newmodel.CMPS(varmet,:);...
                     newmodel.CMPS(cnstmet,:)];
end
if isfield(newmodel,'K')
    newmodel.K = [newmodel.K(varmet,:);...
                  newmodel.K(cnstmet,:)];
elseif ~isempty(newpvec)
    newpvec.K = [newpvec.K(varmet,:);newpvec.K(cnstmet,:)];
end

if isfield(newmodel,'Klb')
    newmodel.Klb = [newmodel.Klb(varmet,:);...
                    newmodel.Klb(cnstmet,:)];
elseif ~isempty(newpvec)
    newpvec.Klb = [newpvec.Klb(varmet,:);...
                   newpvec.Klb(cnstmet,:)];
end

if isfield(newmodel,'Kub')
    newmodel.Kub = [newmodel.Kub(varmet,:);...
                    newmodel.Kub(cnstmet,:)];             
elseif ~isempty(newpvec)
    newpvec.Kub = [newpvec.Kub(varmet,:);...
                   newpvec.Kub(cnstmet,:)];
end

if isfield(newmodel,'KIact')
    newmodel.KIact = [newmodel.KIact(varmet,:);...
                      newmodel.KIact(cnstmet,:)];
elseif ~isempty(newpvec)
    newpvec.KIact = [newpvec.KIact(varmet,:);...
                     newpvec.KIact(cnstmet,:)];
end

if isfield(newmodel,'KIihb')
    newmodel.KIihb = [newmodel.KIihb(varmet,:);...
                      newmodel.KIihb(cnstmet,:)];
elseif ~isempty(newpvec)
    newpvec.KIihb = [newpvec.KIihb(varmet,:);...
                     newpvec.KIihb(cnstmet,:)];
end

if isfield(newpvec,'Vmax')
    newpvec.Vmax = newpvec.Vmax;
end
if isfield(newpvec,'kcat_fwd')
    newpvec.kcat_fwd = newpvec.kcat_fwd;
end
if isfield(newpvec,'kcat_bkw')
    newpvec.kcat_bkw = newpvec.kcat_bkw;
end
if isfield(newpvec,'delGr')
    newpvec.delGr = newpvec.delGr;
end
if isfield(newpvec,'feasible')
    newpvec.feasible = newpvec.feasible;
end

if isfield(newmodel,'enzs')
    newmodel.enzs = newmodel.enzs;
end
if isfield(newmodel,'rxns')
    newmodel.rxns = newmodel.rxns;
end
if isfield(newmodel,'Vss')
    newmodel.Vss = newmodel.Vss;
end
if isfield(newmodel,'delSGr')
    newmodel.delSGr = newmodel.delSGr;
end
if isfield(newmodel,'delGlb')
    newmodel.delGlb = newmodel.delGlb;
end
if isfield(newmodel,'delGub')
    newmodel.delGub = newmodel.delGub;
end
if isfield(newmodel,'Keq')
    newmodel.Keq = newmodel.Keq;
end
if isfield(newmodel,'rev')
    newmodel.rev = newmodel.rev;
end
if isfield(newmodel,'Vact_ind')
    newmodel.Vact_ind = newmodel.Vact_ind;
end
if isfield(newmodel,'Vihb_ind')
    newmodel.Vihb_ind = newmodel.Vihb_ind;
end
if isfield(newmodel,'MClow')
    newmodel.MClow = [newmodel.MClow(varmet,:);...
                      newmodel.MClow(cnstmet,:)];
end
if isfield(newmodel,'MChigh')
    newmodel.MChigh = [newmodel.MChigh(varmet,:);...
                       newmodel.MChigh(cnstmet,:)];
end
if isfield(newmodel,'b')
    newmodel.b = newmodel.b;
end
if isfield(newmodel,'c')
    newmodel.c = newmodel.c;
end
if isfield(newmodel,'vl')
    newmodel.vl = newmodel.vl;
end
if isfield(newmodel,'vu')
    newmodel.vu = newmodel.vu;
end
if isfield(newmodel,'MolWt')
    newmodel.MolWt = [newmodel.MolWt(varmet,:);...
                      newmodel.MolWt(cnstmet,:)];
end
if isfield(newmodel,'Vuptake')
    newmodel.Vuptake = newmodel.Vuptake;
end
if isfield(newmodel,'rxn_add')
    newmodel.rxn_add = newmodel.rxn_add;
end
if isfield(newmodel,'rxn_excep')
    newmodel.rxn_excep = newmodel.rxn_excep;
end

if ~isempty(mc)
    newmc = [newmc(varmet,:);...
             newmc(cnstmet,:)];
end






% function specifyconstantsasparameters