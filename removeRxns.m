function [newmodel,newpvec] = removeRxns(model,pvec,rmid,kpid)
if nargin < 4
    kpid = [];
end

newmodel = model;
rxn_add = model.rxn_add;
rxn_excep = model.rxn_excep;
id = addToVind(model,rxn_add,rxn_excep);

if isempty(rmid)    
    newpvec = pvec;
    return
else
    if any(ismember(rmid,id))
        rmid(ismember(rmid,id)) = [];
    end
    rtid = setdiff(1:model.nt_rxn,rmid);
end

if ~isempty(kpid)
    rtid = kpid;
end

newmodel.rxns = model.rxns(rtid);
newmodel.enzs = model.enzs(rtid);
newmodel.S = model.S(:,rtid);
newmodel.SI = model.SI(:,rtid);
newmodel.CMPS = model.CMPS(:,rtid);
newmodel.Vss = model.Vss(rtid);
newmodel.delSGr = model.delSGr;
newmodel.delGlb = model.delGlb;
newmodel.delGub = model.delGub;
newmodel.Keq = model.Keq;
newmodel.rev = model.rev;
newmodel.c = model.c(rtid);
newmodel.vl = model.vl(rtid);
newmodel.vu = model.vu(rtid);

newmodel.nt_rxn = length(newmodel.rxns);
newmodel.nt_metab = length(newmodel.mets);

%estimate new indices for Vex,
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel,newmodel.nt_rxn,newmodel.S);

newmodel.Vind = Vind;
newmodel.VFex = VFex;
newmodel.Vex = Vex;
newmodel.bmrxn = bmrxn;

newmodel.n_rxn = length(Vind);

newmodel.Vuptake = model.Vuptake(rtid);

newpvec = pvec;
newpvec.K = pvec.K(:,rtid);
% newpvec.KI = pvec.KI(:,rtid);
newpvec.Klb = pvec.Klb(:,rtid);
newpvec.Kub = pvec.Kub(:,rtid);
newpvec.KIact = pvec.KIact(:,rtid);
newpvec.KIihb = pvec.KIihb(:,rtid);
newpvec.Vmax = pvec.Vmax(rtid);
newpvec.kcat_fwd = pvec.kcat_fwd(rtid);
newpvec.kcat_bkw = pvec.kcat_bkw(rtid);
