%function [dXdt,flag,newdata] = ODEmodel(t,Y,data,pmeter)
%**************************************************************************
%Describing the ODE dXdt = data.S*flux(X,p)
%September 2014
%**************************************************************************
function [dXdt,flag,newdata] = ODEmodel(t,Y,model,pmeter)
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,model,pmeter);
bm_ind = model.bmrxn;
Mbio = strcmpi('biomass',model.Metabolites);
% Mbio = model.S(:,bm_ind)>0;

%glcpts
Mglc = strcmpi('glc[e]',model.Metabolites);
Vglc = find(model.S(Mglc,:)<0);
model.Vex = setdiff(model.Vex,Vglc);
model.V2rct = setdiff(model.V2rct,Vglc);

nr = zeros(4,1);
nr(1) = model.nt_metab;
nr(2) = model.nint_metab;%-length(find(Mbio));
nr(3) = model.nt_rxn;
nr(4) = model.n_rxn;
Vind = model.Vind;

% Vic_exind = data.Vic_exind;
% [~,Vuptake] = find(data.S(:,Vic_exind)>0);
[intm_ind,~] = find(model.S(:,model.Vexind)<0);
% Vuptake = Vic_exind(Vuptake);
% Vexcrt = Vic_exind(Vexcrt);
% Mbio = data.S(:,end);
% Mbio_ind = find(Mbio<0);
% nic_exrxn = length(Vic_exind);
% Vxc_exind = data.Vxc_exind;
% nxc_exrxn = length(Vxc_exind);
dXdt = zeros(nr(1),1);%[Metabolites;Biomass]
% dXdt = cons(dXdt,Y);
flux = zeros(nr(3),1);
%% %Biomass flux
% gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
% gr_flux = 0.001;%model.gmax;%0.8;%h-1

%% Fluxes 
%Fixed Uptake Fluxes
flux(model.Vupind) = model.Vuptake;%[20];
%Transporters
flux = ExFlux(model,Y,flux,model.Vex);
flux = ExFlux(model,Y,flux,model.V2rct);
flux(Vglc) = ConvinienceKinetics(model,pmeter,Y,bm_ind,Vglc);
%Intracellular(Cytosolic)
%flux(Vind) = ConvinienceKinetics();
flux(Vind) = ConvinienceKinetics(model,pmeter,Y,bm_ind,Vind);
%Growth Flux
gr_flux = biomass_flux(model,Y,dXdt,flux);
flux(bm_ind) = gr_flux;
fprintf('Growth rate = %3.4f h-1\n',flux(bm_ind));
% flux(model.Vex) = ConvinienceKinetics(model,pmeter,Y,bm_ind,model.Vex);
plotflux_timecourse(flux,t,model);
plotconc_timecourse(Y,t,model);
%Intracellular(Mitochondria) (Yeast)
% flux(VMit) = ConvinienceKinetics(model,pmeter,Y,bm_ind,VMit);
%% %Intracellular Metabolites
%Cytosolic
%ATP, AMP, ADP
% ec = 0.8;
% ATP = strcmpi('atp[c]',model.Metabolites);
% ADP = strcmpi('adp[c]',model.Metabolites);
% AMP = strcmpi('amp[c]',model.Metabolites);
% AdID = [find(ATP),find(ADP),find(AMP)];
% 
% dXdt(AMP) = model.S(AMP,:)*flux - Y(AMP)*gr_flux;
% dXdt(ATP) = model.S(ATP,:)*flux - Y(ATP)*gr_flux;
% dXdt(ADP) = model.S(ADP,:)*flux - Y(ADP)*gr_flux;
% 
% % Y(AMP) = Y(ATP)*(1/ec-1)+Y(ADP)*(1/(2*ec)-1);
% %NAD+, NADH, NADP, NADPH
% NAD = strcmpi('nad[c]',model.Metabolites);
% NADH = strcmpi('nadh[c]',model.Metabolites);
% NADP = strcmpi('nadp[c]',model.Metabolites);
% NADPH = strcmpi('nadph[c]',model.Metabolites);
% NaID = [find(NAD) find(NADH),find(NADP),find(NADPH)];
% 
% dXdt(NAD) = model.S(NAD,:)*flux - Y(NAD)*gr_flux;
% dXdt(NADH) = model.S(NAD,:)*flux - Y(NADH)*gr_flux;
% dXdt(NADP) = model.S(NAD,:)*flux - Y(NADP)*gr_flux;
% dXdt(NADPH) = model.S(NAD,:)*flux - Y(NADPH)*gr_flux;

% mind = setdiff(1:nr(2),[AdID,NaID]);

dXdt(1:nr(2)) = model.S(1:nr(2),:)*flux - Y(1:nr(2))*gr_flux;

%Biomass
dXdt(Mbio) = Y(Mbio)*gr_flux;



%Extracellular Metabolites
% dXdt(nr(2)+1:nr(1)) = D*(data.M*data.S(1:nr(2),:)*flux-Y(nr(2)+1:nr(1)))-...
%                       data.M*data.S(1:nr(2),:)*flux;
glce = strcmpi('glc[e]',model.Metabolites);
lace = strcmpi('lac[e]',model.Metabolites);
dXdt(glce) = model.S(glce,:)*flux;
dXdt(lace) = model.S(lace,:)*flux;



%Biomass 
% dXdt(end) = data.S(end,:)*flux;

if any(Y(Y<0))
    flag = -1;
else
    flag = 0;
end
newdata = [];
end





