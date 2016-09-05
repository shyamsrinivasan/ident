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
    %determine which group idx belongs to
    for id = 1:length(idx)
        if ismember(idx(id),Vind)
            idflux(id,:) = CKinetics(model,pvec,M,idx(id));
        elseif ismember(idx(id),Vex)
            idflux(id,:) = TKinetics(model,pvec,M,idx(id));                
        elseif ismember(idx(id),VFex)
            idflux(id,:) = EKinetics(model,pvec,M,idx(id));
        elseif idx(id)==find(strcmpi(model.rxns,'ATPM'))
            sbid = model.S(:,idx(id))<0;
            sbid(h2o) = 0;
            idflux(idx(id),:) = pvec.Vmax(tatpm).*18.84.*...
                                 M(sbid,:)./pvec.K(sbid,tatpm)./...
                                (1+M(sbid,:)./pvec.K(sbid,tatpm));                
        elseif idx(id)==model.bmrxn
            idflux(id,:) = 0;
        end
        % ETC fluxes - not vectorized yet
        ETCrxn = {'ATPS4r','NADH16','CYTBD','SUCDi','FRD7'};
        if any(strcmpi(ETCrxn,model.rxns{idx(id)}))
            flux = ETCflux(model,pvec,M(:,ic),flux(:,ic));
            idflux(id,ic) = flux(idx(id));
        end
    end                       
end  

if ~isempty(idx)
    flux = idflux;
end

% if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
%     all(mc(logical(model.S(:,model.bmrxn)>0))>0)
%     flux(model.bmrxn) = 0.1;
% else
%     flux(model.bmrxn) = 0;
% end

