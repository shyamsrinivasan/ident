function flux = iflux(model,pvec,M,flux,idx)
[~,nc] = size(M);
if nargin <5
    idx = [];
else
    idflux = zeros(length(idx),nc);
end
if nargin<4
    flux = zeros(model.nt_rxn,nc);
%     flux = cons(flux,M);
end
% if isfield(model,'rxn_add');
%     rxn_add = model.rxn_add;
% else
%     rxn_add = {};
% end
% if isfield(model,'rxn_excep')
%     rxn_excep = model.rxn_excep;
% else
%     rxn_excep = {};
% end

h2o = strcmpi(model.mets,'h2o[c]');

S = model.S;
Vind = model.Vind;
Vex = model.Vex;
% other fixed exchaged fluxes
VFex = model.VFex;

% Vind = addToVind(model,model.Vind,rxn_add,rxn_excep);
% rxn_excep = union(rxn_excep,model.rxns(Vind));
% Vex = addToVind(model,model.Vex,[],rxn_excep);

%carbon uptake fluxes
% [flux,~,vcup] = CarbonKinetics(model,pvec,mc,flux);

%redox reaction fluxes in vrem
% [flux,~,vred] = RedoxKinetics(model,pvec,mc,flux);

tatpm = strcmpi(model.rxns,'ATPM');
if isempty(idx)
    if ~isempty(Vind)
        flux(Vind,:) = CKinetics(model,pvec,M,Vind);
    end
    if ~isempty(Vex)
        flux(Vex,:) = TKinetics(model,pvec,M,Vex);
    end
    if ~isempty(VFex)
        flux(VFex,:) = EKinetics(model,pvec,M,VFex);
    end
    % atp maintanance    
    if any(tatpm)
        sbid = model.S(:,tatpm)<0;
        sbid(h2o) = 0;
        flux(tatpm,:) = pvec.Vmax(tatpm).*18.84.*M(sbid,:)./pvec.K(sbid,tatpm)./...
                      (1+M(sbid,:)./pvec.K(sbid,tatpm));
    end

    % biomass
        %     if mc(strcmpi(model.mets,'atp[c]'))>0 &&...
        %        mc(strcmpi(model.mets,'h2o[c]'))>0
        %         flux(strcmpi(model.rxns,'atpm')) = 8.39;
        %     else
        %         flux(strcmpi(model.rxns,'atpm')) = 0;
        %     end
        %     if all(mc(logical(model.S(:,model.bmrxn)<0))>1e-5)
        %         flux(model.bmrxn) = model.Vss(model.bmrxn);%0.01;
        %     elseif any(mc(logical(model.S(:,model.bmrxn)<0))<1e-5)
        %         flux(model.bmrxn) = 0;
        %     end
    %         flux(model.bmrxn,ic) = BMKinetics(model,pvec,M,model.bmrxn);
    %         flux(model.bmrxn,ic) = model.Vss(model.bmrxn)/3600;
    %         flux(strcmpi('GLCpts',model.rxns),ic) = 10;      
        
else
    % determine which group idx belongs to
    % vectorize wrt idx
    % intracellular reactions
    Vin_idx = idx(ismember(idx,Vind));
    idx = setdiff(idx,Vin_idx);
    if ~isempty(Vin_idx)
        idflux(Vin_idx,:) = CKinetics(model,pvec,M,Vin_idx);
    end
    
    % transport fluxes
    Vex_idx = idx(ismember(idx,Vex));
    idx = setdiff(idx,Vex_idx);
    if ~isempty(Vex_idx)
        idflux(Vex_idx,:) = TKinetics(model,pvec,M,Vex_idx);
    end
    
    % exchange fluxes
    VFex_idx = idx(ismember(idx,VFex));
    idx = setdiff(idx,VFex_idx);    
    if ~isempty(VFex_idx)
        idflux(VFex_idx,:) = TKinetics(model,pvec,M,VFex_idx);
    end
    
    % atp maintenance
    if any(idx==find(strcmpi(model.rxns,'ATPM')))
        ratpm = idx(idx==find(strcmpi(model.rxns,'ATPM')));
        S(h2o,ratpm)=0;      
        idflux(idx==find(strcmpi(model.rxns,'ATPM')),:) =...
        pvec.Vmax(ratpm).*18.84.*...
        M(S(:,ratpm)<0,:)./pvec.K(S(:,ratpm)<0,ratpm)./...
        (1+M(S(:,ratpm)<0,:)./pvec.K(S(:,ratpm)<0,ratpm));   
    end
    
    % biomass reaction
    if any(idx==model.bmrxn)
        idx(idx==model.bmrxn) = 0;
    end
    
    % ETC fluxes - not vectorized yet
%     ETCrxn = {'ATPS4r','NADH16','CYTBD','SUCDi','FRD7'};
%     if any(strcmpi(ETCrxn,model.rxns{idx(id)}))
%         flux = ETCflux(model,pvec,M(:,ic),flux(:,ic));
%         idflux(id,ic) = flux(idx(id));
%     end                   
end 

if ~isempty(idx)
    flux = idflux;
end

% time factor for fluxes - conversion between seconds <-> hours
flux = flux.*3600;

% if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
%     all(mc(logical(model.S(:,model.bmrxn)>0))>0)
%     flux(model.bmrxn) = 0.1;
% else
%     flux(model.bmrxn) = 0;
% end

