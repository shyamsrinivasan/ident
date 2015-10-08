function solvemetLP(model)
%
newmodel.S = model.S;
newmodel.Keq = model.Keq;
newmodel.Vss = model.Vss;
newmodel.rxns = model.rxns;
% Vind = model.Vind;
% Vex = model.Vex;
% VFex = model.VFex;
% 
% vflux = zeros(size(model.S,2),1);
% vflux([Vind Vex VFex']) = 1;

%remove reactions with zero fluxes
newmodel.S(:,abs(newmodel.Vss)<1e-6) = [];
newmodel.Keq(abs(newmodel.Vss)<1e-6) = [];
newmodel.rxns(abs(newmodel.Vss)<1e-6) = [];
newmodel.Vss(abs(newmodel.Vss)<1e-6) = [];
newmodel.mets = model.mets;
% vflux(abs(newmodelVss)<1e-6) = [];
% vflux = find(vflux);

vspl =  [find(strcmpi(newmodel.rxns,'THD2'))...
        find(strcmpi(newmodel.rxns,'NADH16'))...
        find(strcmpi(newmodel.rxns,'ATPS4r'))...        
        find(strcmpi(newmodel.rxns,'CYTBD'))];
    
vmet = [find(strcmpi(newmodel.mets,'h[c]'))...
        find(strcmpi(newmodel.mets,'h[e]'))...
        find(strcmpi(newmodel.mets,'pi[c]'))];

vh2o =  [find(strcmpi(newmodel.mets,'h2o[c]'))...
         find(strcmpi(newmodel.mets,'h2o[e]'))];   
     
% newmodel.S([vmet vh2o],:) = [];   
% newmodel.mets([vmet vh2o],:) = [];
newmodel.S(vmet,setdiff(1:size(newmodel.S,2),vspl)) = 0;    
newmodel.S(vh2o,:) = [];
newmodel.mets(vh2o,:) = [];

[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
Vind = [Vind vspl find(strcmpi(newmodel.rxns,'GLCpts'))]; 
Vex = setdiff(Vex,Vind);    

newmodel.S(:,[Vex VFex' bmrxn find(strcmpi(newmodel.rxns,'ATPM'))]) = [];
newmodel.Keq([Vex VFex' bmrxn find(strcmpi(newmodel.rxns,'ATPM'))]) = [];
newmodel.Vss([Vex VFex' bmrxn find(strcmpi(newmodel.rxns,'ATPM'))]) = [];
newmodel.rxns([Vex VFex' bmrxn find(strcmpi(newmodel.rxns,'ATPM'))]) = [];

% newmodel.S(strcmpi(model.mets,'h[c]'),:) = [];
% newmodel.S(strcmpi(model.mets,'h[e]'),:) = [];
% newmodel.S(strcmpi(model.mets,'h2o[c]'),:) = [];

%remove empty rows in S
newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);

%nonzero reactions for thermodynamic analysis
[nmet,nrxn] = size(newmodel.S);

%Ax <=b 
A = newmodel.S';
A_ub = repmat(sign(newmodel.Vss),1,nmet).*A;
A_lb = repmat(sign(newmodel.Vss),1,nmet).*(-A);
A = [A_ub;A_lb];

b_ub = sign(newmodel.Vss).*log(newmodel.Keq);
b_lb = sign(newmodel.Vss).*(-(log(1e-8)+log(newmodel.Keq)));
b = [b_ub;b_lb];

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-6);
ub = zeros(nmet,1);
ub(ub==0) = log(3e-2);

%concentrations in M = mole/L, Bennett et al., 2009
lb(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% lb(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
lb(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% lb(strcmpi(newmodel.mets,'pep[c]')) = log(1.46e-4);
% lb(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
lb(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
lb(strcmpi(newmodel.mets,'nadh[c]')) = log(5.45e-5);
lb(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
% lb(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
lb(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
lb(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
lb(strcmpi(newmodel.mets,'cit[c]')) = log(1.10e-3);
lb(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% lb(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
lb(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
lb(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
lb(strcmpi(newmodel.mets,'ac[c]')) = log(0.0497);
% lb(strcmpi(newmodel.mets,'co2[c]')) = log(1.75);
% lb(strcmpi(newmodel.mets,'o2[c]')) = log(1.75);

ub(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% ub(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
ub(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% ub(strcmpi(newmodel.mets,'pep[c]')) = log(1.46e-4);
% ub(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
ub(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
% ub(strcmpi(newmodel.mets,'nadh[c]')) = log(5.45e-5);
% ub(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
ub(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
ub(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
ub(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
ub(strcmpi(newmodel.mets,'cit[c]')) = log(1.10e-3);
ub(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% ub(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
ub(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
ub(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
ub(strcmpi(newmodel.mets,'ac[c]')) = log(0.0497);
ub(strcmpi(newmodel.mets,'co2[c]')) = log(3.0);
ub(strcmpi(newmodel.mets,'o2[c]')) = log(1.6);


%try only pyk
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
vrpe = strcmpi(newmodel.rxns,'rpe');
vrpi = strcmpi(newmodel.rxns,'rpi');
vtkt1 = strcmpi(newmodel.rxns,'tkt1');
vtala = strcmpi(newmodel.rxns,'tala');
vtkt2 = strcmpi(newmodel.rxns,'tkt2');
vpdh = strcmpi(newmodel.rxns,'pdh');
vakgd = strcmpi(newmodel.rxns,'akgdh');
vsucs = strcmpi(newmodel.rxns,'sucoas');
vfrd7 = strcmpi(newmodel.rxns,'frd7');
vfum = strcmpi(newmodel.rxns,'fum');
vmdh = strcmpi(newmodel.rxns,'mdh');
vpta = strcmpi(newmodel.rxns,'ptar');
vack = strcmpi(newmodel.rxns,'ackr');
vnadh16 = strcmpi(newmodel.rxns,'nadh16');
vnadtrd = strcmpi(newmodel.rxns,'nadtrhd');
vthd2 = strcmpi(newmodel.rxns,'THD2');
vatps = strcmpi(newmodel.rxns,'ATPS4r');
vppc = strcmpi(newmodel.rxns,'ppc');

list = [find(vpyk) find(vpts) find(vpgi)...
        find(vpfk) find(vfba) find(vtpi)...
        find(vgapd) find(vpgk) find(vpgm)...
        find(veno) find(vrpe) find(vrpi)...
        find(vtkt1) find(vtala) find(vtkt2)...
        find(vpdh) find(vakgd) find(vsucs)...
        find(vfrd7) find(vfum) find(vmdh)...
        find(vpta) find(vack) find(vppc)...
        find(vnadh16) find(vnadtrd) find(vthd2)... 
        ];

A = A_ub(list,:);%;A_lb(vpyk,:)];
b = b_ub(list);%;b_lb(vpyk)];
cprod = sparse(1,nmet);
% cprod = sparse(1,find(strcmpi(newmodel.mets,'atp[c]')),1,1,nmet);

[x,xobj,flag] = cplexlp(-cprod(:),A,b,[],[],lb,ub);


if flag > 0
    %reverse check
    Acheck = repmat(sign(newmodel.Vss(list)),1,nmet).*A;
    bcheck = sign(newmodel.Vss(list)).*b;
    fprintf('We have a solution');
    if A*x<=b
        fprintf('feasible')
    end    
    xub = zeros(1,1);
else
    fprintf('We do not have a solution');
end
