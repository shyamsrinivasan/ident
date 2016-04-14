function [newmodel,newpvec,newmc] = addremoveMets(model,rmvmet,pvec,mc)
if nargin<4
    mc = [];
else
    newmc = mc;
end
if nargin<3
    newpvec = struct([]);
else
    newpvec = pvec;
end
if nargin<2
    rmvmet={};
end

newmodel = model;
if ~isempty(rmvmet)
    list = setdiff(model.mets,rmvmet);
    list = cellfun(@(x)strcmpi(x,model.mets),list,'UniformOutput',false);
    list = cell2mat(cellfun(@(x)find(x),list,'UniformOutput',false));
end

newmodel.S = newmodel.S(list,:);
newmodel.mets = newmodel.mets(list);
newmodel.CMPS = newmodel.CMPS(list,:);
newmodel.SI = newmodel.SI(list,:);
newmodel.MClow = newmodel.MClow(list);
newmodel.MChigh = newmodel.MChigh(list);
newmodel.MolWt = newmodel.MolWt(list);
newmodel.b = newmodel.b(list);

if ~isempty(rmvmet)
    if ~isempty(newpvec)
        newpvec.K = newpvec.K(list,:);
        newpvec.Klb = newpvec.Klb(list,:);
        newpvec.Kub = newpvec.Kub(list,:);
        newpvec.KIact = newpvec.KIact(list,:);
        newpvec.KIihb = newpvec.KIihb(list,:);        
    end
end

% remove empty columns from matrices of nmet x nrxn and vector nrxn x 1 
newmodel.rxns = newmodel.rxns(logical(sum(logical(newmodel.S),1)));
newmodel.SI = newmodel.SI(:,logical(sum(logical(newmodel.S),1)));
newmodel.CMPS = newmodel.CMPS(:,logical(sum(logical(newmodel.S),1)));
newmodel.Keq = newmodel.Keq(logical(sum(logical(newmodel.S),1)));
newmodel.Vss = newmodel.Vss(logical(sum(logical(newmodel.S),1)));
newmodel.delSGr = newmodel.delSGr(logical(sum(logical(newmodel.S),1)));
newmodel.delGlb = newmodel.delGlb(logical(sum(logical(newmodel.S),1)));
newmodel.delGub = newmodel.delGub(logical(sum(logical(newmodel.S),1)));
newmodel.rev = newmodel.rev(logical(sum(logical(newmodel.S),1)));
newmodel.c = newmodel.c(logical(sum(logical(newmodel.S),1)));
newmodel.vl = newmodel.vl(logical(sum(logical(newmodel.S),1)));
newmodel.vu = newmodel.vu(logical(sum(logical(newmodel.S),1)));
newmodel.Vuptake = newmodel.Vuptake(logical(sum(logical(newmodel.S),1)));

% remove corresponding columns in pvec
if ~isempty(newpvec)
    newpvec.K = newpvec.K(:,logical(sum(logical(newmodel.S),1)));
    newpvec.Klb = newpvec.Klb(:,logical(sum(logical(newmodel.S),1)));
    newpvec.Kub = newpvec.Kub(:,logical(sum(logical(newmodel.S),1)));
    newpvec.KIact = newpvec.KIact(:,logical(sum(logical(newmodel.S),1)));
    newpvec.KIihb = newpvec.KIihb(:,logical(sum(logical(newmodel.S),1))); 
    newpvec.Vmax = newpvec.Vmax(logical(sum(logical(newmodel.S),1)));
    newpvec.kcat_fwd = newpvec.kcat_fwd(logical(sum(logical(newmodel.S),1)));
    newpvec.kcat_bkw = newpvec.kcat_bkw(logical(sum(logical(newmodel.S),1)));
    newpvec.delGr = newpvec.delGr(logical(sum(logical(newmodel.S),1)));    
end

newmodel.S = newmodel.S(:,logical(sum(logical(newmodel.S),1)));

% separate cytosolic from external metabolites
[newmodel,newpvec] = separate_cex(newmodel,newpvec);

% calculate new reaction indices
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel,size(newmodel.S,2),newmodel.S);
newmodel.Vind = Vind;
newmodel.VFex = VFex;
newmodel.Vex = Vex;
newmodel.bmrxn = bmrxn;


if ~isempty(mc)
    newmc = newmc(list);
end