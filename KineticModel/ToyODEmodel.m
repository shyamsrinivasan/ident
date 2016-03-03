function [dXdt,flag,newdata] = ToyODEmodel(t,mc,data,model,pvec)

nt_m = model.nt_metab;
nin_m = model.nint_metab;
dXdt = zeros(nt_m,1);

%% Fluxes
flux = iflux(model,pvec,mc.*model.imc);

%% Metabolites
%Cytosolic
dXdt(1:nin_m) = (1./model.imc(1:nin_m)).*(model.S(1:nin_m,:)*flux);

%h[c] is assume constant
hc = strcmpi('h[c]',model.mets);
dXdt(hc) = 0;

%Extracellular
dXdt(nin_m+1:nt_m) = 0;%

%% staus check for CVODE in SUNDIALS TB
if any(mc<0)
    flag = -1;
else
    flag = 0;
end
newdata = data;
newdata.flux = flux;
newdata.Y = mc;
newdata.t = t;