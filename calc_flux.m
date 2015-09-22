function flux = calc_flux(model,pmeter,MC,flux,EC)
%Calculate intial and final flux from concentrations using Conveninece
%Kinetics
%Shyam 2014
useVmax = 0;
if nargin < 5
    useVmax = 1;
end

if nargin < 4 || isempty(flux)
    flux = zeros(model.nt_rxn,1);
end

%transport fluxes 
Vex = model.Vex;
%glucose
Vglc = strcmpi(model.rxns,'glcpts');
if any(Vex == find(Vglc))
    Vex(Vex==find(Vglc)) = [];
    if useVmax
        flux(Vglc) =...
        ConvinienceKinetics(model,pmeter,MC,find(Vglc));
    else
        flux(Vglc) = ExFlux(model,pmeter,MC,flux,find(Vglc),[]);
    end
end
%oxygen
Vo2 = strcmpi(model.rxns,'O2t');
if any(Vex == find(Vo2))
    Vex(Vex==find(Vo2)) = [];
    if useVmax
        flux = ExFlux(model,pmeter,MC,flux,find(Vo2),[]);
    end
end
%other exchange fluxes
Vnadh = find(strcmpi(model.rxns,'NADH16'));
if any(Vex==Vnadh)
    Vex(Vex==Vnadh) = [];
    if useVmax
        flux(Vnadh) =...
        ConvinienceKinetics(model,pmeter,MC,Vnadh);
    end
end

Vthd2 = find(strcmpi(model.rxns,'THD2')); 
if any(Vex==Vthd2)
    Vex(Vex==Vthd2) = [];
    if useVmax
        flux(Vthd2) =...
        ConvinienceKinetics(model,pmeter,MC,Vthd2);
    end
end

Vatps4r = find(strcmpi(model.rxns,'ATPS4r')); 
if any(Vex==Vatps4r)
    Vex(Vex==Vatps4r) = [];
    if useVmax
        flux(Vatps4r) =...
        ConvinienceKinetics(model,pmeter,MC,Vatps4r);
    end
end

Vcyt = find(strcmpi(model.rxns,'CYTBD'));
if any(Vex==Vcyt)
    Vex(Vex==Vcyt) = [];
    if useVmax
        flux(Vcyt) =...
        ConvinienceKinetics(model,pmeter,MC,Vcyt);
    end
end

%other fluxes
Vind = model.Vind;%intracellular fluxes
VFex = model.VFex;%exchange fluxes
Vuptake = model.Vuptake;%fixed uptake fluxes for exchange reactions
if useVmax
    %Transport fluxes
    flux = ExFlux(model,pmeter,MC,flux,Vex,[]);    
    %Intracellular Fluxes
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,Vind,[]);
%     flux(Vex) = pmeter.Vmax(Vex).*MC(int_ind);    
else
    %Transporters
%     if ~isempty(Vex) && ~isempty(int_ind)
%         flux(Vex) = 1000*EC(Vex).*MC(int_ind);
%     end
    flux = ExFlux(model,MC,flux,Vex,[],EC);       
    %Intracellular fluxes
    flux(Vind) = ConvinienceKinetics(model,pmeter,MC,Vind,EC);    
end
%exchange fluxes
% Voth = setdiff(1:model.nt_rxn,VFex);
for ivfex = 1:length(VFex)
%     vexh2o = strcmpi(model.rxns,'exH2O');
%     if VFex(ivfex) == find(vexh2o)   
%         %consuming h2o
%         Vch2o = model.S(strcmpi(model.mets,'h2o[e]'),Voth)<0;
%         Vph2o = model.S(strcmpi(model.mets,'h2o[e]'),Voth)>0;
%         Vuptake(VFex(ivfex)) = sum(flux(Voth(Vch2o)))+...
%                                sum(flux(Voth(Vph2o)));
%     end
%     vexpi = strcmpi(model.rxns,'exPI');
%     if VFex(ivfex) == find(vexpi)
%         %reaction consuming pi
%         Vconpi = model.S(strcmpi(model.mets,'pi[e]'),Voth)<0;
%         Vpropi = model.S(strcmpi(model.mets,'pi[e]'),Voth)>0;
%         Vuptake(VFex(ivfex)) = sum(flux(Voth(Vconpi)))+...
%                                sum(flux(Voth(Vpropi)));
%     end
%     vexh = strcmpi(model.rxns,'exH');
%     if VFex(ivfex) == find(vexh)
%         %consuming h
%         Vconh = model.S(strcmpi(model.mets,'h[e]'),Voth)<0;
%         Vproh = model.S(strcmpi(model.mets,'h[e]'),Voth)>0;
%         Vuptake(VFex(ivfex)) = sum(flux(Voth(Vconh)))+...
%                                sum(flux(Voth(Vproh)));
%     end
%     vexco2 = strcmpi(model.rxns,'exCO2');
%     if VFex(ivfex) == find(vexco2)
%         %consuming h
%         Vconco = model.S(strcmpi(model.mets,'co2[e]'),Voth)<0;
%         Vproco = model.S(strcmpi(model.mets,'co2[e]'),Voth)>0;
%         Vuptake(VFex(ivfex)) = sum(flux(Voth(Vconco)))+...
%                                sum(flux(Voth(Vproco)));
%     end
%     flux(VFex(ivfex)) = f(met[e]) - Vuptake(VFex(ivfex));    
    flux(VFex(ivfex)) = -Vuptake(VFex(ivfex));
end
%biomass Flux
bmind = model.bmrxn;
if isfield(model,'gmax');
    flux(bmind) = model.gmax;
else
    flux(bmind) = biomass_flux(model,MC,[],flux);
end

% if any(Vupind == find(Vglc))
%     Vupind(Vupind==find(Vglc)) = [];
%     if useVmax
%         flux(Vglc) =...
%         ConvinienceKinetics(model,pmeter,MC,bmind,find(Vglc));
%     else
%     end
% else    
%     flux = ExFlux(model,MC,flux,Vupind,'mm');        
% end

% 
% flux(VFup) = flux(Vup);
% flux(VFex) = flux(Vdn);

%biomass exchange
% Vbmex = model.S(strcmpi('biomass',model.mets),:)<0;
% flux(Vbmex) = flux(bmind);

% % gr_flux = 0.8*prod(Y(Mbio_ind)./([.8;.1]+Y(Mbio_ind)));
% gr_flux = model.gmax;%0.8;%h-1
% % gr_flux = biomass_flux(model,Y,dXdt,flux);
% % if isfield(data,'Y') && isfield(data,'t') && isfield(data,'flux')
% %     gr_flux = biomass(model,Y,data.Y,t,data.t,data.flux,flux);
% % else
% %     gr_flux = biomass(model,Y,[],t,[],[],flux);
% % end

% flux(VFup) = flux(Vup);
% flux(VFex) = flux(Vdn);
return
