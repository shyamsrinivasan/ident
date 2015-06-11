%function [flux,vcontr] = ConvinienceKinetics(model,pmeter,MC,bm_ind)
%**************************************************************************
%September 2014
%Changed parameters to sparse matrices
%January 2015 
%**************************************************************************
function [flux,vcontr,varargout] = ConvinienceKinetics(model,pmeter,MC,bm_rxn,Vind,EC)

if nargin < 6
    use_Vmax = 1;    
else
    use_Vmax = 0;
end
%Convinience Kinetics for intracellular fluxes in Metabolism
% n_rxn = length(Vind);
[~,n_rxn] = size(model.S);
% nt_rxn = model.nt_rxn;
% nt_metab = model.nt_metab;
% nint_metab = model.nint_metab;

vthermo_ = zeros(n_rxn,1);
vsat_ = zeros(n_rxn,1);
vact_ = zeros(n_rxn,1);
vihb_ = zeros(n_rxn,1);
vreg_ = zeros(n_rxn,1);
gamma = zeros(n_rxn,1);
flux = zeros(n_rxn,1);
% Vind = model.Vind;

% vcontr = zeros(nt_rxn,1);
jmetab = 0;
kact = 0;
kihb = 0;
h2oid = strcmpi(model.Metabolites,'h2o[c]');
%Non-Exchange Reactions
for irxn = 1:length(Vind)
    subsind = model.S(:,Vind(irxn)) < 0;%substrate
    prodind = model.S(:,Vind(irxn)) > 0;
    %Remove h2o from consideration of kinetics - Non limiting
    if any(subsind(h2oid)) || any(prodind(h2oid))
        subsind(h2oid) = 0;
        prodind(h2oid) = 0;
    end
    if any(subsind) || any(prodind)
        nsubs = length(find(subsind));
        nprod = length(find(prodind));
        if isempty(bm_rxn)||(~isempty(bm_rxn) && Vind(irxn) ~= bm_rxn) 
            %not biomass reaction
            subs = MC(subsind);
            prud = MC(prodind);
            s_subs = -(model.S(subsind, Vind(irxn)));
            s_prod = model.S(prodind, Vind(irxn));
            Ksubs = pmeter.K(subsind,Vind(irxn));
            Kprod = pmeter.K(prodind,Vind(irxn));      
            %Thermodynamics
            if model.Keq(Vind(irxn)) ~= 0%avoid divide by zero error
                if prod(subs.^s_subs) ~= 0%avoid divide by zero error
                    vthermo_(Vind(irxn)) =...
                    1 - prod(prud.^s_prod)/...
                    (prod(subs.^s_subs)*model.Keq(Vind(irxn)));
                    gamma(Vind(irxn)) = prod(subs.^s_subs)/...
                    (prod(prud.^s_prod)*model.Keq(Vind(irxn)));
                else
                    vthermo_(Vind(irxn)) = 0;
                end
            end
            %Saturation
            numrsat = prod(subs.^s_subs);
            drsatsubs = prod((1 + subs./Ksubs).^s_subs);        
            drsatprod = prod((1 + prud./Kprod).^s_prod);        
            vsat_(Vind(irxn)) = numrsat/(drsatsubs*drsatprod - 1);
            jmetab = jmetab + nsubs + nprod;    
        end
    else
        vsat_(Vind(irxn)) = 1;
        vthermo_(Vind(irxn)) = 1;
    end
end
%Regulatory Contribution
%Activators
Vact_ind = model.Vact_ind;
if ~isempty(model.Vact_ind)
    Vact = setdiff(model.Vact_ind,Vind);
    if isempty(Vact)        
        nact_rxn = length(Vact_ind);
        for iact = 1:nact_rxn    
            act_mind = model.SI(:,Vact_ind(iact))>0;
            nact = length(find(act_mind));
            actmc = MC(act_mind);
            s_act = model.SI(act_mind,Vact_ind(iact));
            KIact = pmeter.KIact(act_mind,Vact_ind(iact));
        %     KIact = pmeter.KIact(kact+1:kact+nact);
            vact_(Vact_ind(iact)) = prod(((actmc./KIact).^s_act)./(1 + (actmc./KIact).^s_act));
            kact = kact + nact;    
        end        
    end
    vact_(setdiff(1:n_rxn,Vact_ind)) = 1;
else
    vact_(1:n_rxn) = 1;
end

%Inhibitors
Vinb_ind = model.Vihb_ind;
if ~isempty(model.Vihb_ind)
    Vihb = setdiff(model.Vihb_ind,Vind);
    if isempty(Vihb)        
        nihb_rxn = length(Vinb_ind);
        for jhib = 1:nihb_rxn
            ihb_mind = model.SI(:,Vinb_ind(jhib))<0;
            nihb = length(find(ihb_mind));
            ihbmc = MC(ihb_mind);    
            s_inhib = -(model.SI(ihb_mind,Vinb_ind(jhib)));%elements in SI are -ve 
            KIinb = pmeter.KIihb(ihb_mind,Vinb_ind(jhib));
        %     KIinb = pmeter.KIihb(kihb+1:kihb+nihb); 
            vihb_(Vinb_ind(jhib)) = prod(1./(1 + (ihbmc./KIinb).^s_inhib));  
            kihb = kihb + nihb; 
        end        
    end
    vihb_(setdiff(1:n_rxn,Vinb_ind)) = 1;
else
    vihb_(1:n_rxn) = 1;
end
vreg_(1:n_rxn) = vact_.*vihb_;            

%Net Flux - Intracellular
% model.Kcat(1:n_rxn) = 3000;%s-1
if use_Vmax
%     flux(Vind) = pmeter.Vmax(Vind).*vthermo_(Vind).*vsat_(Vind).*vreg_(Vind);
    flux(Vind) = pmeter.Vmax(Vind).*vthermo_(Vind).*vsat_(Vind).*vreg_(Vind);
else
    flux(Vind) = model.Kcat(Vind).*EC(Vind).*vthermo_(Vind).*vsat_(Vind).*vreg_(Vind);
end
flux = flux(Vind);
%Exchange Fluxes are determined outside in the ODE based on
%1. Material balance
%2. Kinetics
%3. Thermodynamics
vcontr = vthermo_.*vsat_.*vreg_;
varargout{1} = vthermo_;
varargout{2} = vsat_;
varargout{3} = vreg_;
% [mindx,rindx] = find(model.S);
% Jtherm = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
% Jsat = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
% relJtherm = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);
% relJsat = sparse(rindx,mindx,ones(length(rindx),1),nrxn,nmetab);

% [regind,regrind] = find(model.SI);
% Jreg = sparse(regrind,regind,ones(length(regrind),1),nrxn,nmetab);
% relJreg = sparse(regrind,regind,ones(length(regrind),1),nrxn,nmetab);

%Notes: numrsat, drsat, vact and vinhib are declared as arrays just for debugging. They
%will evetually be double scalars that are calculated for each reaction

for irxn = 1:n_rxn
% Initialization
    %indices for each reaction j
%     subsind = model.S(:,irxn) < 0;%substrate
%     prodind = model.S(:,irxn) > 0;%product
%     actind = model.SI(:,irxn) > 0;%activator
%     inhibind = model.SI(:,irxn) < 0;%inhibitor
%     
%     nsubs = length(find(subsind));
%     nprod = length(find(prodind));
%     nact = length(find(actind));
%     ninb = length(find(inhibind));  
%     
%     if ~isempty(biomassind)
%         if irxn ~= find(biomassind)
%             %Substrates & Products

%             %Activators/Inhibitors
    
  
    %Calculate Jacobian (called for each reaction)
    %Inputs: Ksubs, Kprod, KIact, KIinb, s_subs, s_prod, s_act, s_inhib,
    %        subs, pruds
    %Outputs: Jthermo, Jsat, Jreg
    %Jacobian = Jthermo + Jsat + Jreg
    
      
%     [Jsatrx,Jthermrx,Jregrx,...
%      relJsatrx,relJthermrx,relJregrx] = calcJacobian(subs,pruds,act,inhib,...
%                                           Ksubs,Kprod,KIact,KIinb,...
%                                           s_subs,s_prod,s_act,s_inhib,...
%                                           numrsat(irxn),drsatsubs,drsatprod,...
%                                           vact(irxn),vinhib(irxn),...
%                                           gamma(irxn));                          



end
% J = Jsat + Jtherm + Jreg;   


end