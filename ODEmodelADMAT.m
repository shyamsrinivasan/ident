function dXdt = ODEmodelADMAT(mc,model)
%function dXdt = ODEmodelADMAT(mc,model,pvec)
%**************************************************************************
%Describing the ODE dXdt = data.S*flux(X,p)
%September 2014
%**************************************************************************
% callODEmodel = @(t,Y,data)ODEmodel(t,Y,model,pmeter);
bmind = model.bmrxn;
pvec = model.pvec;
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

dXdt = zeros(nt_m,1);%[Metabolites;Biomass]
dXdt = cons(dXdt,mc);
% if isfield(data,'flux')
%     flux = data.flux;
% else

% end

%% Fluxes 
% imc = cons(model.imc,mc);
flux = ifluxADMAT(model,pvec,mc.*model.imc);

% %% %Biomass flux - growth rate
% % gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
% gr_flux = model.gmax;%0.8;%h-1
% gr_flux = biomass_flux(model,mc,[],flux);
% plot(t,gr_flux,'LineStyle','none','Marker','o');
% hold on
% % if isfield(data,'Y') && isfield(data,'t') && isfield(data,'flux')
% %     gr_flux = biomass(model,Y,data.Y,t,data.t,data.flux,flux);
% % else
% %     gr_flux = biomass(model,Y,[],t,[],[],flux);
% % end
% mu = flux(bmind);

%plot time course concentrations and flux during integration
% plotflux_timecourse(flux,t,model,[nad16 atps cyt pit act]);


%% %Intracellular Metabolites
%Cytosolic
dXdt(1:nin_m) = (1./model.imc(1:nin_m)).*(model.S(1:nin_m,:)*flux);%-mu*mc(1:nin_m);%
% dXdt(1:nin_m) = -mu*mc(1:nin_m);%
dXdt(nin_m+1:nt_m) = zeros(length(nin_m+1:nt_m),1);
% dXdt(nin_m+1:nt_m) = (1./model.imc(nin_m+1:nt_m)).*(model.S(nin_m+1:nt_m,:)*zeros(model.nt_rxn,1));
%model.S(nin_m+1:nt_m,:)*flux;%-mu*Y(nin_m+1:nt_m);

% plotconc_timecourse(dXdt,t,model,[hc he pic pie]);
idx = find(mc<0);
% if t > 1e-7
%     dbstop in ODEmodel.m at 136
% end
% if any(idx)
%     dbstop in ODEmodel.m at 145
%     fprintf('%d %3.6g %d\n',length(idx),t,idx(:));
%     
% %     return
% end

% 
% %Biomass
% dXdt(Mbio) = model.S(Mbio,:)*flux*Y(Mbio);

%Extracellular Metabolites
% dXdt(strcmpi('A[e]',model.mets)) = model.S(strcmpi('A[e]',model.mets),:)*flux;
% dXdt(nr(2)+1:nr(1)) = D*(data.M*data.S(1:nr(2),:)*flux-Y(nr(2)+1:nr(1)))-...
%                       data.M*data.S(1:nr(2),:)*flux;

%Biomass 
% dXdt(end) = data.S(end,:)*flux;
% 
% newdata = data;
% newdata.flux = flux;
% newdata.Y = mc;
% newdata.t = t;
end





