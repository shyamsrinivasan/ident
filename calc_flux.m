function flux = calc_flux(model,pmeter,MC,flux,EC)
%Calculate intial and final flux from concentrations using Conveninece
%Kinetics
%Shyam 2014
useVmax = 0;
if nargin < 5
    useVmax = 1;
end
nt_rxn = model.nt_rxn;
if nargin < 4
    flux = zeros(nt_rxn,1);
end

Vind = model.Vind;
Vupind = model.Vupind;
Vex = model.Vex;
Vup = model.Vup;
Vdn = model.Vdn;
VFup = model.VFup;
VFex = model.VFex;
bmind = model.bmrxn;

if ~isempty(Vex)
    [int_ind,~] = find(model.S(:,Vex)<0);
else
    int_ind = [];
end

%Uptake Fluxes
flux(Vupind) = 20;%ExFlux(model,MC,flux,Vupind,'mm');

if useVmax
    %Transporters
    flux = ExFlux(model,MC,flux,Vex,[]);
    %Intracellular Fluxes
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,bmind,Vind);
%     flux(Vex) = pmeter.Vmax(Vex).*MC(int_ind);    
else
    %Transporters
%     if ~isempty(Vex) && ~isempty(int_ind)
%         flux(Vex) = 1000*EC(Vex).*MC(int_ind);
%     end
    flux = ExFlux(model,MC,flux,Vex,[],EC);   
    %Intracellular fluxes
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,bmind,Vind,EC);     
end

flux(VFup) = flux(Vup);
flux(VFex) = flux(Vdn);
%Biomass Flux
if isfield(model,'gmax');
    flux(bmind) = model.gmax;
else
    flux(bmind) = biomass_flux(model,MC,[],flux);
end
% % gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
% gr_flux = model.gmax;%0.8;%h-1
% % gr_flux = biomass_flux(model,Y,dXdt,flux);
% % if isfield(data,'Y') && isfield(data,'t') && isfield(data,'flux')
% %     gr_flux = biomass(model,Y,data.Y,t,data.t,data.flux,flux);
% % else
% %     gr_flux = biomass(model,Y,[],t,[],[],flux);
% % end

flux(VFup) = flux(Vup);
flux(VFex) = flux(Vdn);
return
