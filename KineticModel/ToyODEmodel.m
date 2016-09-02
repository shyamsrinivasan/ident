function [dXdt,flag,newdata] = ToyODEmodel(t,mc,data,model,pvec)

nvar = length(mc);
dXdt = zeros(nvar,1);
Pimc = model.Pimc;

%% Fluxes
flux = iflux(model,pvec,[mc.*model.imc;Pimc]);

%% Metabolites
% Cytosolic
dXdt(1:nvar) = (1./model.imc).*(model.S(1:nvar,:)*(flux*3600));

% h[c] is assume constant
% hc = strcmpi('h[c]',model.mets);
% dXdt(hc) = 0;

% change pi[c]
% pic = strcmpi('pi[c]',model.mets);
% dXdt(pic) = dXdt(pic)+0.001;

%Extracellular
% dXdt(nin_m+1:nt_m) = 0;%

%% staus check for CVODE in SUNDIALS TB
% if any(mc<0)
%     flag = -1;
% else
    flag = 0;
% end
newdata = data;
newdata.flux = flux;
newdata.Y = mc;
newdata.t = t;