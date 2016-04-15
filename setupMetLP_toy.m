function newmodel = setupMetLP_toy(model)

%create new model structure
newmodel.S = model.S;
newmodel.Keq = model.Keq;
newmodel.Vss = model.Vss;
newmodel.rxns = model.rxns;
newmodel.CMPS = model.CMPS;

%remove reactions with zero fluxes
newmodel.S(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Keq(abs(newmodel.Vss)<1e-7) = [];
newmodel.rxns(abs(newmodel.Vss)<1e-7) = [];
newmodel.CMPS(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Vss(abs(newmodel.Vss)<1e-7) = [];
newmodel.mets = model.mets;

%special reactions - reactions where h[c], h[e],pi[c] affect 
%thermodynamic equilibrium 
vspl = [];

%metabolites that do not affect thermodynamic equilibrium   
vmet= [];

%find all exchnage and transport reactions in newmodel
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
% Vind = ToColumnVector(Vind);
% VFex = ToColumnVector(VFex);
% Vex = ToColumnVector(Vex);
Vind = [Vind vspl find(strcmpi(newmodel.rxns,'Ain'))]; 
Vex = setdiff(Vex,Vind);     

%remove h[c], h[e],pi[c]
%remove from all reactions except vspl
if ~isempty(vmet)
    newmodel.S(vmet,:) = [];
    newmodel.mets(vmet) = [];
end

%remove h2o[c]
vh2o = [find(strcmpi(newmodel.mets,'h2o[c]'))...
        find(strcmpi(newmodel.mets,'h2o[e]'))];  
    
vhe = [];

newmodel.S([vh2o vhe],:) = [];
newmodel.mets([vhe vh2o]) = []; 

newmodel.S(:,[Vex VFex bmrxn]) = [];
newmodel.Keq([Vex VFex bmrxn]) = [];
newmodel.Vss([Vex VFex bmrxn]) = [];
newmodel.rxns([Vex VFex bmrxn]) = [];

%remove empty rows in S
newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);

vatpm = find(strcmpi(newmodel.rxns,'ATPM'));  

newmodel.Vss = newmodel.Vss(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.rxns = newmodel.rxns(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.Keq = newmodel.Keq(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.S = newmodel.S(:,setdiff(1:size(newmodel.S,2),vatpm));

%nonzero reactions for thermodynamic analysis
nmet = size(newmodel.S,1);

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-8);
ub = zeros(nmet,1);
ub(ub==0) = log(3e-2);

%concentrations in M = mole/L, Bennett et al., 2009
for imet = 1:length(newmodel.mets)
    tfm = strcmpi(newmodel.mets{imet},model.mets);
    if any(tfm)
        lb(imet) = log(model.MClow(tfm));
        ub(imet) = log(model.MChigh(tfm));
    end
end
lb(strcmpi(newmodel.mets,'A[e]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'A[c]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'B[c]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'C[c]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'D[c]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'E[c]')) = log(0.2);
% 
ub(strcmpi(newmodel.mets,'A[e]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'A[c]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'B[c]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'C[c]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'D[c]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'E[c]')) = log(0.2);

knwn_id = zeros(nmet,1);
for imet = 1:nmet
    if lb(imet) == ub(imet)
        knwn_id(imet)=1;
    end
end

%setup original problem
%Ax <=b 
A = newmodel.S';
A_ub = A(:,~logical(knwn_id));
newmodel.A = A_ub;
newmodel.A_kn = A(:,logical(knwn_id));

b_ub = log(newmodel.Keq)-A(:,logical(knwn_id))*lb(logical(knwn_id));
newmodel.b = b_ub;

newmodel.x = lb(logical(knwn_id));%known concentrations
newmodel.lb = lb(~logical(knwn_id));
newmodel.ub = ub(~logical(knwn_id));
newmodel.mets_kn = newmodel.mets(logical(knwn_id));
newmodel.mets = newmodel.mets(~logical(knwn_id));










