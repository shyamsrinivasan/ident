function [LPmax,LPmin] = solvemetLP(newmodel,prxnid)
if nargin<2
    prxnid = 0;
end

%temporary call
% newmodel = setupMetLP(model);
nmet = size(newmodel.S,1);

%try only pyk
% vpyk = strcmpi(newmodel.rxns,'pyk');
% vpts = strcmpi(newmodel.rxns,'glcpts');
% vpgi = strcmpi(newmodel.rxns,'pgi');
% vpfk = strcmpi(newmodel.rxns,'pfk');
% vfba = strcmpi(newmodel.rxns,'fba');
% vtpi = strcmpi(newmodel.rxns,'tpi');
% vgapd = strcmpi(newmodel.rxns,'gapd');
% vpgk = strcmpi(newmodel.rxns,'pgk');
% vpgm = strcmpi(newmodel.rxns,'pgm');
% veno = strcmpi(newmodel.rxns,'eno');
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
% vnadh16 = strcmpi(newmodel.rxns,'nadh16');
% vnadtrd = strcmpi(newmodel.rxns,'nadtrhd');
% vthd2 = strcmpi(newmodel.rxns,'THD2');
% vatps = strcmpi(newmodel.rxns,'ATPS4r');
% vcytb = strcmpi(newmodel.rxns,'CYTBD');
% vppc = strcmpi(newmodel.rxns,'ppc');

% list = [find(vpyk) find(vpts) find(vpgi)...
%         find(vpfk) find(vfba) find(vtpi)...
%         find(vgapd) find(vpgk) find(vpgm)...
%         find(veno) find(vrpe) find(vrpi)...
%         find(vtkt1) find(vtala) find(vtkt2)...
%         find(vpdh) find(vakgd) find(vsucs)...
%         find(vfrd7) find(vfum) find(vmdh)...
%         find(vpta) find(vack) find(vppc)...
%         find(vnadh16) find(vnadtrd) find(vthd2)...         
%         find(vatps) find(vcytb)];

A = newmodel.A;
b = newmodel.b;
lb = newmodel.lb;
ub = newmodel.ub;

if prxnid
    cprod = sparse(1,prxnid,1,1,nmet);
else
    cprod = sparse(1,nmet);
end

%maximization
[x,xobj,flag] = cplexlp(-cprod(:),A,b,[],[],lb,ub);
if flag>0    
    LPmax.x = x;
    LPmax.obj = xobj;
end
LPmax.flag = flag;

%minimization
if prxnid
    [x,xobj,flag] = cplexlp(cprod(:),A,b,[],[],lb,ub);
    if flag>0    
        LPmin.x = x;
        LPmin.obj = xobj;
    end
    LPmin.flag = flag;
else
    LPmin = struct([]);
end




    
