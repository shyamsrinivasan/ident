function [newmodel,newpvec,newmc] = separate_cex(model,pvec,mc)
if nargin<3
    mc = [];
    newmc = mc;
end
if nargin<2
    pvec = [];
    newpvec = pvec;
end
% Separate External & Internal mets
% exter_mind = ~cellfun('isempty',regexp(model.mets,'\w(?:xt)$'));
exter_mind = ~cellfun('isempty',regexp(model.mets,'\w(?:\[e\])$'));
inter_mind = ~exter_mind;
bm_ind = strcmpi(model.mets,'Biomass[c]');
if any(bm_ind)
    model.mets{bm_ind} = 'Biomass';
end
inter_mind(bm_ind)=0;

newmodel.mets = [model.mets(inter_mind,1);...
                 model.mets(exter_mind,1);...
                 model.mets(bm_ind,1)];   
newmodel.S = [model.S(inter_mind,:);model.S(exter_mind,:);model.S(bm_ind,:)];
newmodel.SI = [model.SI(inter_mind,:);model.SI(exter_mind,:);model.SI(bm_ind,:)];
newmodel.next_metab = length(find(exter_mind));
newmodel.nint_metab = length(find(inter_mind));
if isfield(model,'CMPS')
    newmodel.CMPS = [model.CMPS(inter_mind,:);...
                     model.CMPS(exter_mind,:);...
                     model.CMPS(bm_ind,:)];
end
if isfield(model,'K')
    newmodel.K = [model.K(inter_mind,:);...
                  model.K(exter_mind,:);...
                  model.K(bm_ind,:)];
elseif ~isempty(pvec)
    newpvec.K = [pvec.K(inter_mind,:);pvec.K(exter_mind,:);pvec.K(bm_ind,:)];
end

if isfield(model,'Klb')
    newmodel.Klb = [model.Klb(inter_mind,:);...
                    model.Klb(exter_mind,:);...
                    model.Klb(bm_ind,:)];
elseif ~isempty(pvec)
    newpvec.Klb = [pvec.Klb(inter_mind,:);...
                   pvec.Klb(exter_mind,:);...
                   pvec.Klb(bm_ind,:)];
end

if isfield(model,'Kub')
    newmodel.Kub = [model.Kub(inter_mind,:);...
                    model.Kub(exter_mind,:);...
                    model.Kub(bm_ind,:)];             
elseif ~isempty(pvec)
    newpvec.Kub = [pvec.Kub(inter_mind,:);...
                   pvec.Kub(exter_mind,:);...
                   pvec.Kub(bm_ind,:)];
end

if isfield(model,'KIact')
    newmodel.KIact = [model.KIact(inter_mind,:);...
                      model.KIact(exter_mind,:);...
                      model.KIact(bm_ind,:)];
elseif ~isempty(pvec)
    newpvec.KIact = [pvec.KIact(inter_mind,:);...
                     pvec.KIact(exter_mind,:);...
                     pvec.KIact(bm_ind,:)];
end

if isfield(model,'KIihb')
    newmodel.KIihb = [model.KIihb(inter_mind,:);...
                      model.KIihb(exter_mind,:);...
                      model.KIihb(bm_ind,:)];
elseif ~isempty(pvec)
    newpvec.KIihb = [pvec.KIihb(inter_mind,:);...
                     pvec.KIihb(exter_mind,:);...
                     pvec.KIihb(bm_ind,:)];
end
if ~isempty(pvec)
    if isfield(pvec,'Vmax')
        newpvec.Vmax = pvec.Vmax;
    end
    if isfield(pvec,'kfwd')
        newpvec.kfwd = pvec.kfwd;
    end
    if isfield(pvec,'krev')
        newpvec.krev = pvec.krev;
    end
    if isfield(pvec,'delGr')
        newpvec.delGr = pvec.delGr;
    end
    if isfield(pvec,'feasible')
        newpvec.feasible = pvec.feasible;
    end
end
if isfield(model,'enzs')
    newmodel.enzs = model.enzs;
end
if isfield(model,'rxns')
    newmodel.rxns = model.rxns;
end
if isfield(model,'Vss')
    newmodel.Vss = model.Vss;
end
if isfield(model,'delSGr')
    newmodel.delSGr = model.delSGr;
end
if isfield(model,'delGlb')
    newmodel.delGlb = model.delGlb;
end
if isfield(model,'delGub')
    newmodel.delGub = model.delGub;
end
if isfield(model,'Keq')
    newmodel.Keq = model.Keq;
end
if isfield(model,'rev')
    newmodel.rev = model.rev;
end
if isfield(model,'Vact_ind')
    newmodel.Vact_ind = model.Vact_ind;
end
if isfield(model,'Vihb_ind')
    newmodel.Vihb_ind = model.Vihb_ind;
end
if isfield(model,'MClow')
    newmodel.MClow = [model.MClow(inter_mind,:);...
                      model.MClow(exter_mind,:);...
                      model.MClow(bm_ind,:)];
end
if isfield(model,'MChigh')
    newmodel.MChigh = [model.MChigh(inter_mind,:);...
                       model.MChigh(exter_mind,:);...
                       model.MChigh(bm_ind,:)];
end
if isfield(model,'b')
    newmodel.b = model.b;
end
if isfield(model,'c')
    newmodel.c = model.c;
end
if isfield(model,'vl')
    newmodel.vl = model.vl;
end
if isfield(model,'vu')
    newmodel.vu = model.vu;
end
if isfield(model,'MolWt')
    newmodel.MolWt = [model.MolWt(inter_mind,:);...
                      model.MolWt(exter_mind,:);...
                      model.MolWt(bm_ind,:)];
end
if isfield(model,'Vuptake')
    newmodel.Vuptake = model.Vuptake;
end
if isfield(model,'rxn_add')
    newmodel.rxn_add = model.rxn_add;
end
if isfield(model,'rxn_excep')
    newmodel.rxn_excep = model.rxn_excep;
end

if ~isempty(mc)
    newmc = [mc(inter_mind,:);...
             mc(exter_mind,:);...
             mc(bm_ind,:)];
end

    
