function [newmodel,newpvec,newmc] = addremoveRxns(model,rmvrxn,pvec,mc)
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
    rmvrxn = {};
end

newmodel = model;

% try only pyk
vpyk = strcmpi(newmodel.rxns,'pyk');
vpts = strcmpi(newmodel.rxns,'glcpts');
vpgi = strcmpi(newmodel.rxns,'pgi');
vpfk = strcmpi(newmodel.rxns,'pfk');
vfba = strcmpi(newmodel.rxns,'fba');
vtpi = strcmpi(newmodel.rxns,'tpi');
vgapd = strcmpi(newmodel.rxns,'gapd');
vpgk = strcmpi(newmodel.rxns,'pgk');
vpgm = strcmpi(newmodel.rxns,'pgm');
veno = strcmpi(newmodel.rxns,'eno');
vpfl = strcmpi(newmodel.rxns,'pfl');
vacald = strcmpi(newmodel.rxns,'acald');
vacld2x = strcmpi(newmodel.rxns,'alcd2x');
% vrpe = strcmpi(newmodel.rxns,'rpe');
% vrpi = strcmpi(newmodel.rxns,'rpi');
% vtkt1 = strcmpi(newmodel.rxns,'tkt1');
% vtala = strcmpi(newmodel.rxns,'tala');
% vtkt2 = strcmpi(newmodel.rxns,'tkt2');
% vpdh = strcmpi(newmodel.rxns,'pdh');
% vakgd = strcmpi(newmodel.rxns,'akgdh');
% vsucs = strcmpi(newmodel.rxns,'sucoas');
% vfrd7 = strcmpi(newmodel.rxns,'frd7');
% vfum = strcmpi(newmodel.rxns,'fum');
% vmdh = strcmpi(newmodel.rxns,'mdh');
% vpta = strcmpi(newmodel.rxns,'ptar');
% vack = strcmpi(newmodel.rxns,'ackr');
vnadh16 = strcmpi(newmodel.rxns,'nadh16');
% vnadtrd = strcmpi(newmodel.rxns,'nadtrhd');
% vthd2 = strcmpi(newmodel.rxns,'THD2');
vatps = strcmpi(newmodel.rxns,'ATPS4r');
vcytb = strcmpi(newmodel.rxns,'CYTBD');
% vppc = strcmpi(newmodel.rxns,'ppc');

if isempty(rmvrxn)
    list = [find(vpyk) find(vpts) find(vpgi)...
            find(vpfk) find(vfba) find(vtpi)...
            find(vgapd) find(vpgk) find(vpgm)...
            find(veno) find(vpfl) find(vacald)...
            find(vacld2x) find(vnadh16) find(vcytb)...
            find(vatps)];% find(vrpe) find(vrpi)...
    %         find(vtkt1) find(vtala) find(vtkt2)...
    %         find(vpdh) find(vakgd) find(vsucs)...
    %         find(vfrd7) find(vfum) find(vmdh)...
    %         find(vpta) find(vack) find(vppc)...
    %         find(vnadh16) find(vnadtrd) find(vthd2)...         
    %         find(vatps) find(vcytb)];
else
    list = setdiff(model.rxns,rmvrxn);
    list = cellfun(@(x)strcmpi(x,model.rxns),list,'UniformOutput',false);
    list = cell2mat(cellfun(@(x)find(x),list,'UniformOutput',false));
end

newmodel.S = newmodel.S(:,list);
newmodel.SI = newmodel.SI(:,list);
newmodel.Keq = newmodel.Keq(list);
newmodel.rxns = newmodel.rxns(list);
newmodel.CMPS = newmodel.CMPS(:,list);
newmodel.Vss = newmodel.Vss(list);
if isfield(newmodel,'delSGr')
    newmodel.delSGr = newmodel.delSGr(list);
end
if isfield(newmodel,'delGlb')    
    newmodel.delGlb = newmodel.delGlb(list);
end
if isfield(newmodel,'delGub')
    newmodel.delGub = newmodel.delGub(list);
end
if isfield(newmodel,'rev')
    newmodel.rev = newmodel.rev(list);
end
if isfield(newmodel,'c')
    newmodel.c = newmodel.c(list);
end
if isfield(newmodel,'vl')
    newmodel.vl = newmodel.vl(list);
end
if isfield(newmodel,'vu')
    newmodel.vu = newmodel.vu(list);
end
if isfield(newmodel,'Vuptake')
    newmodel.Vuptake = newmodel.Vuptake(list);
end
if ~isempty(rmvrxn)
    if ~isempty(newpvec)
        newpvec.K = newpvec.K(:,list);
        newpvec.Klb = newpvec.Klb(:,list);
        newpvec.Kub = newpvec.Kub(:,list);
        newpvec.KIact = newpvec.KIact(:,list);
        newpvec.KIihb = newpvec.KIihb(:,list); 
        newpvec.Vmax = newpvec.Vmax(list);
        newpvec.kfwd = newpvec.kfwd(list);
        newpvec.krev = newpvec.krev(list);
        if isfield(newpvec,'delGr')
            newpvec.delGr = newpvec.delGr(list);   
        end
    end
end

% remove empty rows from matrices nmet x nrxn or vectors nmet x 1
if isfield(newmodel,'mets')
    newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
end
if isfield(newmodel,'remid')
    if ~isempty(newmodel.remid)
        oldmets = newmodel.mets(newmodel.remid);
        newmetid = cellfun(@(x)strcmpi(newmodel.mets,x),oldmets,'UniformOutput',false);
        newmetid = cellfun(@(x)find(x),newmetid,'UniformOutput',false);
        newmodel.remid = cell2mat(newmetid);
    end
end
if isfield(newmodel,'SI')
    newmodel.SI = newmodel.SI(logical(sum(logical(newmodel.S),2)),:);
end
if isfield(newmodel,'CMPS')
    newmodel.CMPS = newmodel.CMPS(logical(sum(logical(newmodel.S),2)),:);
end
if isfield(newmodel,'MClow')
    newmodel.MClow = newmodel.MClow(logical(sum(logical(newmodel.S),2)));
end
if isfield(newmodel,'MChigh')
    newmodel.MChigh = newmodel.MChigh(logical(sum(logical(newmodel.S),2)));
end
if isfield(newmodel,'MolWt')
    newmodel.MolWt = newmodel.MolWt(logical(sum(logical(newmodel.S),2)));
end

% remove corepsonding rows in pvec
if ~isempty(newpvec)
    newpvec.K = newpvec.K(logical(sum(logical(newmodel.S),2)),:);
    newpvec.Klb = newpvec.Klb(logical(sum(logical(newmodel.S),2)),:);
    newpvec.Kub = newpvec.Kub(logical(sum(logical(newmodel.S),2)),:);
    newpvec.KIact = newpvec.KIact(logical(sum(logical(newmodel.S),2)),:);
    newpvec.KIihb = newpvec.KIihb(logical(sum(logical(newmodel.S),2)),:);     
end

if ~isempty(mc)
    newmc = newmc(logical(sum(logical(newmodel.S),2)));
end

newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);
% calculate new reaction indices
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel,size(newmodel.S,2),newmodel.S);
if isfield(newmodel,'Vind')
    newmodel.Vind = Vind;
end
if isfield(newmodel,'VFex')
    newmodel.VFex = VFex;
end
if isfield(newmodel,'Vex')
    newmodel.Vex = Vex;
end
if isfield(newmodel,'bmrxn')
    newmodel.bmrxn = bmrxn;
end
if isfield(newmodel,'nt_metab')
    newmodel.nt_metab = length(newmodel.mets);
end
if isfield(newmodel,'nt_rxn')
    newmodel.nt_rxn = length(newmodel.rxns);
end

    

