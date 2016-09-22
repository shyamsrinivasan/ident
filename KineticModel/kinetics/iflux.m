function [flux,fluxbm] = iflux(model,pvec,M,flux,idx)
[~,nc] = size(M);
if nargin <5
    idx = [];
else
    idflux = zeros(model.nt_rxn,nc);
end
if nargin<4
    flux = zeros(model.nt_rxn,nc);
%     flux = cons(flux,M);
end

fluxbm = zeros(size(M,1),nc);
h2o = strcmpi(model.mets,'h2o[c]');

S = model.S;
Vind = model.Vind;
Vex = model.Vex;
VFex = model.VFex; % other exchage fluxes
if isfield(model,'bmrxn')
    bmrxn = model.bmrxn;
end

%carbon uptake fluxes
% [flux,~,vcup] = CarbonKinetics(model,pvec,mc,flux);

%redox reaction fluxes in vrem
% [flux,~,vred] = RedoxKinetics(model,pvec,mc,flux);

tatpm = strcmpi(model.rxns,'ATPM');
finid = []; % only when idx is non empty

if isempty(idx)
    if ~isempty(Vind)
        % flux in mmole/Lcw/s
        flux(Vind,:) = CKinetics(model,pvec,M,Vind);
    end
    if ~isempty(Vex)
        % flux in mmole/Lcw/s
        flux(Vex,:) = TKinetics(model,pvec,M,Vex);
    end
    if ~isempty(VFex)
        % flux in mmole/Lc/s
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
    fluxbm = BMKinetics(model,pvec,M,bmrxn);
    flux(bmrxn,:) = 1.0;
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
        % flux in mmole/Lcw/s
        idflux(Vin_idx,:) = CKinetics(model,pvec,M,Vin_idx);
        finid = [finid Vin_idx]; 
    end
    
    % transport fluxes
    Vex_idx = idx(ismember(idx,Vex));
    idx = setdiff(idx,Vex_idx);
    if ~isempty(Vex_idx)
        % flux in mmole/Lcw/s
        idflux(Vex_idx,:) = TKinetics(model,pvec,M,Vex_idx);
        finid = [finid Vex_idx]; 
    end
    
    % exchange fluxes
    VFex_idx = idx(ismember(idx,VFex));
    idx = setdiff(idx,VFex_idx);    
    if ~isempty(VFex_idx)
        % flux in mmole/Lc/s
        idflux(VFex_idx,:) = EKinetics(model,pvec,M,VFex_idx);
        finid = [finid VFex_idx]; 
    end
    
    % atp maintenance
    if ~isempty(idx)
        if ~isempty(find(strcmpi(model.rxns,'ATPM'),1))
            if any(idx==find(strcmpi(model.rxns,'ATPM')))
                ratpm = idx(idx==find(strcmpi(model.rxns,'ATPM')));
                S(h2o,ratpm)=0;      
                idflux(idx==find(strcmpi(model.rxns,'ATPM')),:) =...
                pvec.Vmax(ratpm).*18.84.*...
                M(S(:,ratpm)<0,:)./pvec.K(S(:,ratpm)<0,ratpm)./...
                (1+M(S(:,ratpm)<0,:)./pvec.K(S(:,ratpm)<0,ratpm));   
                finid = [finid ratpm]; 
            end
        end
    end
    
    % biomass reaction
    if ~isempty(idx)
        if ~isempty(model.bmrxn)
            if any(idx==model.bmrxn)
                idflux(idx==model.bmrxn) = 0;
                finid = [finid idx(idx==model.bmrxn)]; 
            end
        end
    end
    
    % ETC fluxes - not vectorized yet
%     ETCrxn = {'ATPS4r','NADH16','CYTBD','SUCDi','FRD7'};
%     if any(strcmpi(ETCrxn,model.rxns{idx(id)}))
%         flux = ETCflux(model,pvec,M(:,ic),flux(:,ic));
%         idflux(id,ic) = flux(idx(id));
%     end                   
end 

if ~isempty(finid)
    flux = idflux(finid);
end

% time factor for fluxes - conversion between seconds -> hours
flux = flux.*3600;
fluxbm = fluxbm.*3600;

% if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
%     all(mc(logical(model.S(:,model.bmrxn)>0))>0)
%     flux(model.bmrxn) = 0.1;
% else
%     flux(model.bmrxn) = 0;
% end

