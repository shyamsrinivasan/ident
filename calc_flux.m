function flux = calc_flux(model,pmeter,MC,EC)
%Calculate intial and final flux from concentrations using Conveninece
%Kinetics
%Shyam 2014
useVmax = 0;
if nargin < 4
    useVmax = 1;
end

nt_rxn = model.nt_rxn;
Vind = model.Vind;
Vupind = model.Vupind;
Vexind = model.Vexind;
if ~isempty(Vexind)
    [int_ind,~] = find(model.S(:,Vexind)<0);
else
    int_ind = [];
end
Vex = model.Vex;
bm_ind = model.bmrxn;
flux = zeros(nt_rxn,1);

flux(Vupind) = model.Vuptake;
if useVmax
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,bm_ind,Vind);
    flux(Vexind) = pmeter.Vmax(Vexind).*MC(int_ind);
    flux = ExFlux(model,MC,flux,Vex);
else
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,bm_ind,Vind,EC);
    if ~isempty(Vexind) && ~isempty(int_ind)
        flux(Vexind) = 1000*EC(Vexind).*MC(int_ind);
    end
    flux = ExFlux(model,MC,flux,Vex,EC);    
end
flux(bm_ind) = model.gmax;
return
