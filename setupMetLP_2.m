function newmodel = setupMetLP_2(model)
%determines concentrations for reactions with delG = 0 and Vss = 0

%create new model structure
newmodel.S = [model.S(:,abs(model.Vss)>1e-7) model.S(:,abs(model.Vss)<1e-7)];
newmodel.Keq = [model.Keq(abs(model.Vss)>1e-7);model.Keq(abs(model.Vss)<1e-7)];
newmodel.rxns = [model.rxns(abs(model.Vss)>1e-7);model.rxns(abs(model.Vss)<1e-7)];
newmodel.Vss = [model.Vss(abs(model.Vss)>1e-7);model.Vss(abs(model.Vss)<1e-7)];
newmodel.mets = model.mets;


%special reactions - reactions where h[c], h[e],pi[c] affect 
%thermodynamic equilibrium 
vspl =  [find(strcmpi(newmodel.rxns,'THD2'))...
        find(strcmpi(newmodel.rxns,'NADH16'))...
        find(strcmpi(newmodel.rxns,'ATPS4r'))... 
        find(strcmpi(newmodel.rxns,'CYTBD'))];

%metabolites that do not affect thermodynamic equilibrium   
vmet = [find(strcmpi(newmodel.mets,'h[c]'))...
        find(strcmpi(newmodel.mets,'h[e]'))...
        find(strcmpi(newmodel.mets,'pi[c]'))];
    
%remove h2o[c]
vh2o =  [find(strcmpi(newmodel.mets,'h2o[c]'))...
         find(strcmpi(newmodel.mets,'h2o[e]'))]; 
    
%find all exchnage and transport reactions in newmodel
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
Vind = ToColumnVector(Vind);
VFex = ToColumnVector(VFex);
Vex = ToColumnVector(Vex);
Vind = [Vind vspl find(strcmpi(newmodel.rxns,'GLCpts'))]; 
Vex = setdiff(Vex,Vind);   

%remove h[c], h[e],pi[c]
newmodel.S([vmet vh2o],:) = [];
newmodel.mets([vmet vh2o]) = [];

newmodel.S(:,[Vex VFex bmrxn]) = [];
newmodel.Keq([Vex VFex bmrxn]) = [];
newmodel.Vss([Vex VFex bmrxn]) = [];
newmodel.rxns([Vex VFex bmrxn]) = [];     

%remove empty rows in S
newmodel.mets = newmodel.mets(logical(sum(logical(newmodel.S),2)));
newmodel.S = newmodel.S(logical(sum(logical(newmodel.S),2)),:);

%zero reactions only
nmet = size(newmodel.S,1);

%bounds
lb = zeros(nmet,1);
lb(lb==0) = log(1e-7);
ub = zeros(nmet,1);
ub(ub==0) = log(3e-2);

%concentrations in M = mole/L, Bennett et al., 2009
lb(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% lb(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
% lb(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% lb(strcmpi(newmodel.mets,'pep[c]')) = log(1.46e-4);
% lb(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
lb(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
lb(strcmpi(newmodel.mets,'nadh[c]')) = log(5.45e-5);
lb(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
% lb(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
lb(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
lb(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
% lb(strcmpi(newmodel.mets,'cit[c]')) = log(1.10e-3);
lb(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% lb(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
lb(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
lb(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
lb(strcmpi(newmodel.mets,'ac[c]')) = log(0.0497);
lb(strcmpi(newmodel.mets,'co2[c]')) = log(1.0);
% lb(strcmpi(newmodel.mets,'o2[c]')) = log(1.75);
% lb(strcmpi(newmodel.mets,'h[c]')) = log(1e-7);
lb(strcmpi(newmodel.mets,'h2o[c]')) = log(55.0);

ub(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% ub(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
% ub(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% ub(strcmpi(newmodel.mets,'pep[c]')) = log(1.46e-4);
% ub(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
ub(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
% ub(strcmpi(newmodel.mets,'nadh[c]')) = log(5.45e-5);
% ub(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
ub(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
ub(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
ub(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
% ub(strcmpi(newmodel.mets,'cit[c]')) = log(1.10e-3);
% ub(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% ub(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
% ub(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
% ub(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
ub(strcmpi(newmodel.mets,'ac[c]')) = log(0.0497);
ub(strcmpi(newmodel.mets,'co2[c]')) = log(3.0);
ub(strcmpi(newmodel.mets,'o2[c]')) = log(1.6);
% ub(strcmpi(newmodel.mets,'h[c]')) = log(1e-7);
ub(strcmpi(newmodel.mets,'h2o[c]')) = log(55.0);


Aeq = xmodel.S';
beq = log(xmodel.Keq);

A = newmodel.S';
A_ub = repmat(sign(newmodel.Vss),1,nmet).*A;

b_ub = sign(newmodel.Vss).*log(newmodel.Keq);

cprod = sparse(1,nmet);

%maximization
[x,xobj,flag] = cplexlp(-cprod(:),A_ub,b_ub,Aeq,beq,lb,ub);
if flag>0
    LPmax.x=x;
    LPmax.obj=xobj;
end
LPmax.flag = flag;

