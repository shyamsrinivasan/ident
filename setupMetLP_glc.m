function newmodel = setupMetLP_glc(model,mc)
if nargin<2
    mc = [];
end
%% %create new model structure
newmodel.S = model.S;
newmodel.Keq = model.Keq;
newmodel.Vss = model.Vss;
newmodel.rxns = model.rxns;
newmodel.CMPS = model.CMPS;

%% %remove reactions with zero fluxes
newmodel.S(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Keq(abs(newmodel.Vss)<1e-7) = [];
newmodel.rxns(abs(newmodel.Vss)<1e-7) = [];
newmodel.CMPS(:,abs(newmodel.Vss)<1e-7) = [];
newmodel.Vss(abs(newmodel.Vss)<1e-7) = [];
newmodel.mets = model.mets;
if ~isempty(mc)
    newmodel.mc = mc;
end

%% %special reactions - reactions where h[c], h[e],pi[c] affect 
%thermodynamic equilibrium 
vspl =  [find(strcmpi(newmodel.rxns,'CYTBD'))];

%% %metabolites that do not affect thermodynamic equilibrium   
vmet= [];

%% %find all exchnage and transport reactions in newmodel
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
Vind = [Vind vspl find(strcmpi(newmodel.rxns,'GLCpts'))]; 
Vex = setdiff(Vex,Vind);     

%% %remove h[c], h[e],pi[c]
%remove from all reactions except vspl
if ~isempty(vmet)
    newmodel.S(vmet,:) = [];
    newmodel.mets(vmet) = [];
    if ~isempty(mc)
        newmodel.mc(vmet) = [];
    end
    % [~,cl] = find(newmodel.CMPS(vmet,:));
%     vspl = unique([vspl]);% ToColumnVector(cl)]);
%     for irxn = 1:length(Vind)
%         if ~ismember(Vind(irxn),vspl)
%             if any(newmodel.S(vmet(logical(newmodel.S(vmet,Vind(irxn)))),Vind(irxn)))
%                 newmodel.S(vmet(logical(newmodel.S(vmet,Vind(irxn)))),Vind(irxn)) = 0;
%             end
%         end
%     end
end

%% %remove h2o[c]
vh2o = [find(strcmpi(newmodel.mets,'h2o[c]'))...
        find(strcmpi(newmodel.mets,'h2o[e]'))];  

vhe = find(strcmpi(newmodel.mets,'h[e]'));    
vhc = find(strcmpi(newmodel.mets,'h[c]'));    
     
% newmodel.S(vh2o,setdiff(1:size(newmodel.S,2),vatpm)) = 0;      
newmodel.S([vh2o vhe vhc],:) = [];
newmodel.mets([vhe vhc vh2o]) = []; 
if ~isempty(mc)
    newmodel.mc([vhe vhc vh2o]) = []; 
end

newmodel.S(:,[Vex VFex bmrxn]) = [];
newmodel.Keq([Vex VFex bmrxn]) = [];
newmodel.Vss([Vex VFex bmrxn]) = [];
newmodel.rxns([Vex VFex bmrxn]) = [];

%% %remove empty rows in S
newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);
if ~isempty(mc)
    newmodel.mc = newmodel.mc(logical(sum(logical(newmodel.S),2)));
end

vatpm = find(strcmpi(newmodel.rxns,'ATPM'));  

newmodel.Vss = newmodel.Vss(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.rxns = newmodel.rxns(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.Keq = newmodel.Keq(setdiff(1:size(newmodel.S,2),vatpm));
newmodel.S = newmodel.S(:,setdiff(1:size(newmodel.S,2),vatpm));

%% %nonzero reactions for thermodynamic analysis
nmet = size(newmodel.S,1);

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-8);
ub = zeros(nmet,1);
ub(ub==0) = log(3e-2);

%concentrations in M = mole/L, Bennett et al., 2009
lb(strcmpi(newmodel.mets,'o2[c]')) = log(mc(strcmpi(model.mets,'o2[c]')));
lb(strcmpi(newmodel.mets,'q8[c]')) = log(mc(strcmpi(model.mets,'q8[c]')));
lb(strcmpi(newmodel.mets,'q8h2[c]')) = log(mc(strcmpi(model.mets,'q8h2[c]')));
lb(strcmpi(newmodel.mets,'mal[c]')) = log(mc(strcmpi(model.mets,'mal[c]')));
lb(strcmpi(newmodel.mets,'fum[c]')) = log(mc(strcmpi(model.mets,'fum[c]')));
lb(strcmpi(newmodel.mets,'nadph[c]')) = log(mc(strcmpi(model.mets,'nadph[c]')));
lb(strcmpi(newmodel.mets,'nadp[c]')) = log(mc(strcmpi(model.mets,'nadp[c]')));
lb(strcmpi(newmodel.mets,'succ[c]')) = log(mc(strcmpi(model.mets,'succ[c]')));
lb(strcmpi(newmodel.mets,'pyr[c]')) = log(mc(strcmpi(model.mets,'pyr[c]')));
lb(strcmpi(newmodel.mets,'co2[c]')) = log(mc(strcmpi(model.mets,'co2[c]')));

ub(strcmpi(newmodel.mets,'o2[c]')) = log(mc(strcmpi(model.mets,'o2[c]')));
ub(strcmpi(newmodel.mets,'q8[c]')) = log(mc(strcmpi(model.mets,'q8[c]')));
ub(strcmpi(newmodel.mets,'q8h2[c]')) = log(mc(strcmpi(model.mets,'q8h2[c]')));
ub(strcmpi(newmodel.mets,'mal[c]')) = log(mc(strcmpi(model.mets,'mal[c]')));
ub(strcmpi(newmodel.mets,'fum[c]')) = log(mc(strcmpi(model.mets,'fum[c]')));
ub(strcmpi(newmodel.mets,'nadph[c]')) = log(mc(strcmpi(model.mets,'nadph[c]')));
ub(strcmpi(newmodel.mets,'nadp[c]')) = log(mc(strcmpi(model.mets,'nadp[c]')));
ub(strcmpi(newmodel.mets,'succ[c]')) = log(mc(strcmpi(model.mets,'succ[c]')));
ub(strcmpi(newmodel.mets,'pyr[c]')) = log(mc(strcmpi(model.mets,'pyr[c]')));
ub(strcmpi(newmodel.mets,'co2[c]')) = log(mc(strcmpi(model.mets,'co2[c]')));


knwn_id = zeros(nmet,1);
for imet = 1:nmet
    if lb(imet) == ub(imet)
        knwn_id(imet)=1;
    end
end

%% %setup original problem
%Ax <=b 
A = newmodel.S';
A_ub = A(:,~logical(knwn_id));
newmodel.A = A_ub;
newmodel.A_kn = A(:,logical(knwn_id));

b_ub = log(newmodel.Keq)-A(:,logical(knwn_id))*lb(logical(knwn_id));
newmodel.b = b_ub;

if ~isempty(mc)
    newmodel.mc = mc(~logical(knwn_id));
end
newmodel.x_kn = lb(logical(knwn_id));%known concentrations
newmodel.lb = lb(~logical(knwn_id));
newmodel.ub = ub(~logical(knwn_id));
newmodel.mets_kn = newmodel.mets(logical(knwn_id));
newmodel.mets = newmodel.mets(~logical(knwn_id));



