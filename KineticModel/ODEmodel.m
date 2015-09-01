%function [dXdt,flag,newdata] = ODEmodel(t,Y,data,pmeter)
%**************************************************************************
%Describing the ODE dXdt = data.S*flux(X,p)
%September 2014
%**************************************************************************
function [dXdt,flag,newdata] = ODEmodel(t,Y,data,model,pmeter)
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,model,pmeter);
bmind = model.bmrxn;
% Mbio = strcmpi('biomass',model.mets);
% % Mbio = model.S(:,bm_ind)>0;
% nt_rxn = model.nt_rxn;
% nin_rxn = model.n_rxn;

% Vind = model.Vind;
% Vex = model.Vex;
% VFex = model.VFex;

nt_m = model.nt_metab;
nin_m = model.nint_metab;
nex_m = model.next_metab;


% nr = zeros(4,1);
% nr(1) = model.nt_metab;
% nr(2) = model.nint_metab;%-length(find(Mbio));
% nr(3) = model.nt_rxn;
% nr(4) = model.n_rxn;


% Vic_exind = data.Vic_exind;
% [~,Vuptake] = find(data.S(:,Vic_exind)>0);
% [intm_ind,~] = find(model.S(:,model.Vexind)<0);
% Vuptake = Vic_exind(Vuptake);
% Vexcrt = Vic_exind(Vexcrt);
% Mbio = data.S(:,end);
% Mbio_ind = find(Mbio<0);
% nic_exrxn = length(Vic_exind);
% Vxc_exind = data.Vxc_exind;
% nxc_exrxn = length(Vxc_exind);
dXdt = zeros(nt_m,1);%[Metabolites;Biomass]
% dXdt = cons(dXdt,Y);
% if isfield(data,'flux')
%     flux = data.flux;
% else

% end

%% Fluxes 
%initial flux
flux = calc_flux(model,pmeter,Y);

% %% %Biomass flux - growth rate
% % gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
% gr_flux = model.gmax;%0.8;%h-1
% % gr_flux = biomass_flux(model,Y,dXdt,flux);
% % if isfield(data,'Y') && isfield(data,'t') && isfield(data,'flux')
% %     gr_flux = biomass(model,Y,data.Y,t,data.t,data.flux,flux);
% % else
% %     gr_flux = biomass(model,Y,[],t,[],[],flux);
% % end
mu = flux(bmind);

%plot time course concentrations and flux during integration
% plotflux_timecourse(flux,t,model)
% plotconc_timecourse(Y,t,model)
%% %Intracellular Metabolites
%Cytosolic
dXdt(1:nin_m) = model.S(1:nin_m,:)*flux-mu*Y(1:nin_m);
dXdt(nin_m+1:nt_m) = model.S(nin_m+1:nt_m,:)*flux;%-mu*Y(nin_m+1:nt_m);

%ATP, AMP, ADP
% ec = 0.8;
% ATP = strcmpi('atp[c]',model.mets);
% ADP = strcmpi('adp[c]',model.mets);
% AMP = strcmpi('amp[c]',model.mets);
% AdID = [find(ATP),find(ADP),find(AMP)];
% 
% dXdt(AMP) = 0;
% dXdt(ATP) = model.S(ATP,:)*flux - Y(ATP)*mu;
% dXdt(ADP) = model.S(ADP,:)*flux - Y(ADP)*mu;
% 
% Y(AMP) = Y(ATP)*(1/ec-1)+Y(ADP)*(1/(2*ec)-1);
% %NAD+, NADH, NADP, NADPH
% NAD = strcmpi('nad[c]',model.mets);
% NADH = strcmpi('nadh[c]',model.mets);
% NADP = strcmpi('nadp[c]',model.mets);
% NADPH = strcmpi('nadph[c]',model.mets);
% NaID = [find(NAD) find(NADH),find(NADP),find(NADPH)];
% 
% dXdt(NAD) = model.S(NAD,:)*flux - Y(NAD)*mu;
% dXdt(NADH) = model.S(NAD,:)*flux - Y(NADH)*mu;
% dXdt(NADP) = model.S(NAD,:)*flux - Y(NADP)*mu;
% dXdt(NADPH) = model.S(NAD,:)*flux - Y(NADPH)*mu;
% 
% mind = setdiff(1:nr(2),[AdID,NaID]);
% mind = setdiff(mind,find(Mbio));
% 
% dXdt(mind) = model.S(mind,:)*flux - Y(mind)*mu;
% 
% %Biomass
% dXdt(Mbio) = model.S(Mbio,:)*flux*Y(Mbio);



%Extracellular Metabolites
% dXdt(strcmpi('A[e]',model.mets)) = model.S(strcmpi('A[e]',model.mets),:)*flux;
% dXdt(nr(2)+1:nr(1)) = D*(data.M*data.S(1:nr(2),:)*flux-Y(nr(2)+1:nr(1)))-...
%                       data.M*data.S(1:nr(2),:)*flux;

%Biomass 
% dXdt(end) = data.S(end,:)*flux;

if any(Y(Y<0))
    flag = -1;
else
    flag = 0;
end
newdata = data;
newdata.flux = flux;
newdata.Y = Y;
newdata.t = t;
end





