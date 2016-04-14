function newmodel = addremoveRxns(model,rmvrxn)
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
newmodel.Keq = newmodel.Keq(list);
newmodel.rxns = newmodel.rxns(list);
newmodel.CMPS = newmodel.CMPS(:,list);
newmodel.Vss = newmodel.Vss(list);


