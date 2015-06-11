%function [dXdt,flag,newdata] = ODEmodel(t,Y,data,pmeter)
%**************************************************************************
%Describing the ODE dXdt = data.S*flux(X,p)
%September 2014
%**************************************************************************
function [dXdt] = ODEmodel_ADMAT(Y,Extra)
data = Extra.model;
pmeter = Extra.pmeter;
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,model,pmeter);
bm_ind = data.bm_ind;
nr = zeros(4,1);
nr(1) = data.nt_metab;
nr(2) = data.nint_metab;
nr(3) = data.nt_rxn;
nr(4) = data.n_rxn;
Vind = data.Vind;
Vic_exind = data.Vic_exind;
Mbio = data.S(:,end);
Mbio_ind = find(Mbio<0);
nic_exrxn = length(Vic_exind);
Vxc_exind = data.Vxc_exind;
nxc_exrxn = length(Vxc_exind);
dXdt = zeros(nr(1)+1,1);%[Metabolites;Biomass]
flux = zeros(nr(3),1);
%Dilution Rate
D = 0.3;
%Net Flux - Exchange(Extracellular)
flux(Vxc_exind) = [-25;25];%fixed exchange flux (uptake and e
%flux(Vxc_exind) = -flux(Vic_exind);
% flux(Vxc_exind) = data.M*data.S(1:nr(2),:)*flux;

%Exchange(Intracellular)
% ex_flux = (data.S(1:nr(2),Vic_exind)*ones(nic_exrxn,1)).*...
%           (data.S(1:nr(2),Vind)*flux(Vind));
% flux(Vic_exind) = ex_flux(1:nic_exrxn);
flux(Vic_exind) = -flux(Vxc_exind);

%Intracellular
%flux(Vind) = ConvinienceKinetics();
flux(Vind) = ConvinienceKinetics(data,pmeter,Y,bm_ind);

%Biomass flux = mumax*(S1/Ks1+S1*S2/Ks2+S2)
gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
flux(end) = (gr_flux-D)*Y(end);

%Intracellular Metabolites
dXdt(1:nr(2)) = data.S(1:nr(2),:)*flux - Y(1:nr(2))*gr_flux;

%Extracellular Metabolites
dXdt(nr(2)+1:nr(1)) = D*(data.M*data.S(1:nr(2),:)*flux-Y(nr(2)+1:nr(1)))-data.M*data.S(1:nr(2),:)*flux;

%Biomass 
dXdt(end) = data.S(end,:)*flux;

flag = 0;
newdata = [];
end





