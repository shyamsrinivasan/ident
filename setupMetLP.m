function newmodel = setupMetLP(model)

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
vspl =  [find(strcmpi(newmodel.rxns,'THD2'))...
        find(strcmpi(newmodel.rxns,'NADH16'))...
        find(strcmpi(newmodel.rxns,'ATPS4r'))... 
        find(strcmpi(newmodel.rxns,'CYTBD'))];

%metabolites that do not affect thermodynamic equilibrium   
vmet= [];%find(strcmpi(newmodel.mets,'h[e]'));%[];
% vmet = [find(strcmpi(newmodel.mets,'h[e]'))...
%         find(strcmpi(newmodel.mets,'h[c]'))...
%         find(strcmpi(newmodel.mets,'pi[c]'))];
     
%find all exchnage and transport reactions in newmodel
[Vind,VFex,Vex,bmrxn] = fluxIndex(newmodel);
Vind = ToColumnVector(Vind);
VFex = ToColumnVector(VFex);
Vex = ToColumnVector(Vex);
Vind = [Vind vspl find(strcmpi(newmodel.rxns,'GLCpts'))]; 
Vex = setdiff(Vex,Vind);     

%remove h[c], h[e],pi[c]
%remove from all reactions except vspl
if ~isempty(vmet)
    newmodel.S(vmet,:) = [];
    newmodel.mets(vmet) = [];
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

%remove h2o[c]
vh2o = [find(strcmpi(newmodel.mets,'h2o[c]'))...
        find(strcmpi(newmodel.mets,'h2o[e]'))];  

vhe = [];%find(strcmpi(newmodel.mets,'h[e]'));    
     
% newmodel.S(vh2o,setdiff(1:size(newmodel.S,2),vatpm)) = 0;      
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
lb(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% lb(strcmpi(newmodel.mets,'fdp[c]')) = log(1e-7);
% lb(strcmpi(newmodel.mets,'2pg[c]')) = log(1e-5);
% lb(strcmpi(newmodel.mets,'3pg[c]')) = log(0.5e-5);
% lb(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% lb(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
% lb(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% lb(strcmpi(newmodel.mets,'pep[c]')) = log(1e-7);
% lb(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
% lb(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
% lb(strcmpi(newmodel.mets,'nadh[c]')) = log(5.45e-5);
% lb(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
% lb(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
% lb(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
% lb(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
lb(strcmpi(newmodel.mets,'cit[c]')) = log(1.10e-3);
lb(strcmpi(newmodel.mets,'icit[c]')) = log(1.0e-5);
% lb(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% lb(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
% lb(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
% lb(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
% lb(strcmpi(newmodel.mets,'ac[c]')) = log(0.0497);
% lb(strcmpi(newmodel.mets,'co2[c]')) = log(1.0);
% lb(strcmpi(newmodel.mets,'o2[c]')) = log(1.75);
% lb(strcmpi(newmodel.mets,'h[c]')) = log(1e-7);
% lb(strcmpi(newmodel.mets,'Xu5p-D[c]')) = log(1e-4);
% lb(strcmpi(newmodel.mets,'ru5p-D[c]')) = log(1e-7);
lb(strcmpi(newmodel.mets,'h2o[c]')) = log(55.0);
% lb(strcmpi(newmodel.mets,'pi[c]')) = log(1e-9);
% lb(strcmpi(newmodel.mets,'h[c]')) = log(1e-8);

ub(strcmpi(newmodel.mets,'glc[e]')) = log(0.2);
% ub(strcmpi(newmodel.mets,'fdp[c]')) = log(1.5e-7);
% ub(strcmpi(newmodel.mets,'2pg[c]')) = log(1e-3);
% ub(strcmpi(newmodel.mets,'f6p[c]')) = log(5e-3);
% ub(strcmpi(newmodel.mets,'atp[c]')) = log(8.13e-3);
% ub(strcmpi(newmodel.mets,'adp[c]')) = log(4.37e-4);
% ub(strcmpi(newmodel.mets,'amp[c]')) = log(2.32e-4);
% ub(strcmpi(newmodel.mets,'pep[c]')) = log(1.46e-4);
% ub(strcmpi(newmodel.mets,'dhap[c]')) = log(3.44e-4);
% ub(strcmpi(newmodel.mets,'nad[c]')) = log(2.32e-3);
% ub(strcmpi(newmodel.mets,'nadh[c]')) = log(2.32e-3);
% ub(strcmpi(newmodel.mets,'nadp[c]')) = log(1.4e-7);
% ub(strcmpi(newmodel.mets,'nadph[c]')) = log(1.1e-4);
% ub(strcmpi(newmodel.mets,'accoa[c]')) = log(5.29e-4);
% ub(strcmpi(newmodel.mets,'coa[c]')) = log(8.8e-5);
ub(strcmpi(newmodel.mets,'cit[c]')) = log(0.25);
% ub(strcmpi(newmodel.mets,'icit[c]')) = log(1.10e-6);
% ub(strcmpi(newmodel.mets,'mal[c]')) = log(1.66e-3);
% ub(strcmpi(newmodel.mets,'fum[c]')) = log(3e-6);
% ub(strcmpi(newmodel.mets,'succ[c]')) = log(3.41e-4);
% ub(strcmpi(newmodel.mets,'succoa[c]')) = log(1.42e-4);
% ub(strcmpi(newmodel.mets,'ac[c]')) = log(1e-3);
ub(strcmpi(newmodel.mets,'co2[c]')) = log(3.0);
ub(strcmpi(newmodel.mets,'o2[c]')) = log(1.6);
% ub(strcmpi(newmodel.mets,'h[c]')) = log(1e-7);
% ub(strcmpi(newmodel.mets,'h2o[c]')) = log(55.0);
% ub(strcmpi(newmodel.mets,'pi[c]')) = log(5e-2);
ub(strcmpi(newmodel.mets,'h[c]')) = log(1e-4);
% ub(strcmpi(newmodel.mets,'h[e]')) = log(2.5e-5);

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
% A_ub = repmat(sign(newmodel.Vss),1,nmet).*A;
newmodel.A = A_ub;
newmodel.A_kn = A(:,logical(knwn_id));

% b_ub = sign(newmodel.Vss).*log(newmodel.Keq);
b_ub = log(newmodel.Keq)-A(:,logical(knwn_id))*lb(logical(knwn_id));
% b_lb = sign(newmodel.Vss).*(-(log(1e-8)+log(newmodel.Keq)));
newmodel.b = b_ub;%b_lb];

newmodel.x = lb(logical(knwn_id));%known concentrations
newmodel.lb = lb(~logical(knwn_id));
newmodel.ub = ub(~logical(knwn_id));
newmodel.mets_kn = newmodel.mets(logical(knwn_id));
newmodel.mets = newmodel.mets(~logical(knwn_id));



% %setup slack problem
% newmodel = setupSlackVariables(newmodel);

% nconstr = length((newmodel.A(:,1)));
% %add slack variables to all constraints
% A_slack = sparse(1:nconstr,1:nconstr,1,nconstr,nconstr);
% A_slack = repmat(sign(newmodel.Vss),1,nconstr).*A_slack;
% 
% newmodel.A = [newmodel.A A_slack];
% 
% lb_slack = zeros(size(newmodel.A,1),1);
% % lb_slack(lb_slack==0) = -Inf;
% ub_slack = zeros(size(newmodel.A,1),1);
% ub_slack(ub_slack==0) = Inf;
% 
% newmodel.lb = [newmodel.lb;lb_slack];
% newmodel.ub = [newmodel.ub;ub_slack];

%add row to A for sum of all slacks
% nvar = size(newmodel.A,2);
% nmets = length(newmodel.mets);
% A_row = sparse(1,nmets+1:nvar,-1,1,nvar);
% newmodel.A = [newmodel.A;A_row];
% 
% %add rows to b
% newmodel.b = [newmodel.b;Inf];
% 
% %add new column for sum of all slacks
% nconstr = size(newmodel.A,1);
% A_col = sparse(nconstr,1,1,nconstr,1);
% newmodel.A = [newmodel.A A_col];
% 
% %add rows to lb and ub
% newmodel.lb = [newmodel.lb;0];
% newmodel.ub = [newmodel.ub;Inf];



