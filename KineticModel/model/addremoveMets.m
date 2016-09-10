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
if isfield(newmodel,'MClow')
    newmodel.MClow = newmodel.MClow(list);
end
if isfield(newmodel,'MChigh')
    newmodel.MChigh = newmodel.MChigh(list);
end
if isfield(newmodel,'MolWt')
    newmodel.MolWt = newmodel.MolWt(list);
end
if isfield(newmodel,'b')
    newmodel.b = newmodel.b(list);
end
if ~isempty(mc)
    newmc = newmc(list);
end

if ~isempty(rmvmet)
    if ~isempty(newpvec)
        if isfield(newpvec,'K')
            newpvec.K = newpvec.K(list,:);
        end
        if isfield(newpvec,'Klb')
            newpvec.Klb = newpvec.Klb(list,:);
        end
        if isfield(newpvec,'Kub')
            newpvec.Kub = newpvec.Kub(list,:);
        end
        if isfield(newpvec,'KIact')
            newpvec.KIact = newpvec.KIact(list,:);
        end
        if isfield(newpvec,'KIihb')
            newpvec.KIihb = newpvec.KIihb(list,:);
        end
    end
end

% remove empty columns from matrices of nmet x nrxn and vector nrxn x 1 
newmodel.rxns = newmodel.rxns(logical(sum(logical(newmodel.S),1)));
newmodel.SI = newmodel.SI(:,logical(sum(logical(newmodel.S),1)));
newmodel.CMPS = newmodel.CMPS(:,logical(sum(logical(newmodel.S),1)));
newmodel.Keq = newmodel.Keq(logical(sum(logical(newmodel.S),1)));
newmodel.Vss = newmodel.Vss(logical(sum(logical(newmodel.S),1)));
if isfield(newmodel,'delSGr')
    newmodel.delSGr = newmodel.delSGr(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'delGlb')
    newmodel.delGlb = newmodel.delGlb(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'delGub')
    newmodel.delGub = newmodel.delGub(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'rev')
    newmodel.rev = newmodel.rev(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'c')
    newmodel.c = newmodel.c(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'vl')
    newmodel.vl = newmodel.vl(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'vu')
    newmodel.vu = newmodel.vu(logical(sum(logical(newmodel.S),1)));
end
if isfield(newmodel,'Vuptake');
    newmodel.Vuptake = newmodel.Vuptake(logical(sum(logical(newmodel.S),1)));
end

% remove corresponding columns in pvec
if ~isempty(newpvec)
    if isfield(newpvec,'K')
        newpvec.K = newpvec.K(:,logical(sum(logical(newmodel.S),1)));
    end
    if isfield(newpvec,'Klb')
        newpvec.Klb = newpvec.Klb(:,logical(sum(logical(newmodel.S),1)));
    end
    if isfield(newpvec,'Kub')
        newpvec.Kub = newpvec.Kub(:,logical(sum(logical(newmodel.S),1)));
    end
    if isfield(newpvec,'KIact')
        newpvec.KIact = newpvec.KIact(:,logical(sum(logical(newmodel.S),1)));
    end
    if isfield(newpvec,'KIihb')
        newpvec.KIihb = newpvec.KIihb(:,logical(sum(logical(newmodel.S),1))); 
    end
    newpvec.Vmax = newpvec.Vmax(logical(sum(logical(newmodel.S),1)));
    newpvec.kfwd = newpvec.kfwd(logical(sum(logical(newmodel.S),1)));
    newpvec.krev = newpvec.krev(logical(sum(logical(newmodel.S),1)));
    newpvec.delGr = newpvec.delGr(logical(sum(logical(newmodel.S),1)));    
end

newmodel.S = newmodel.S(:,logical(sum(logical(newmodel.S),1)));

% separate cytosolic from external metabolites
[newmodel,newpvec,newmc] = separate_cex(newmodel,newpvec,newmc);

% calculate new reaction indices
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel,size(newmodel.S,2),newmodel.S);
newmodel.Vind = Vind;
newmodel.VFex = VFex;
newmodel.Vex = Vex;
newmodel.bmrxn = bmrxn;
newmodel.nt_metab = length(newmodel.mets);
newmodel.nt_rxn = length(newmodel.rxns);


