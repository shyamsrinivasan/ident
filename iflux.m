function flux = iflux(model,pvec,mc,flux,idx)
[~,nc] = size(mc);
if nargin <5
    idx = [];
else
    idflux = zeros(length(idx),nc);
end
if nargin<4
    flux = zeros(model.nt_rxn,nc);
end
if isfield(model,'rxn_add');
    rxn_add = model.rxn_add;
else
    rxn_add = {};
end
if isfield(model,'rxn_excep')
    rxn_excep = model.rxn_excep;
else
    rxn_excep = {};
end

h2o = strcmpi(model.mets,'h2o[c]');

Vind = addToVind(model,model.Vind,rxn_add,rxn_excep);
rxn_excep = union(rxn_excep,model.rxns(Vind));
Vex = addToVind(model,model.Vex,[],rxn_excep);

%carbon uptake fluxes
% [flux,~,vcup] = CarbonKinetics(model,pvec,mc,flux);

%redox reaction fluxes in vrem
% [flux,~,vred] = RedoxKinetics(model,pvec,mc,flux);

%other fixed exchaged fluxes
VFex = model.VFex;

for ic = 1:nc
    if isempty(idx)

        % cytosolic fluxes
        flux(Vind,ic) = CKinetics(model,pvec,mc(:,ic),Vind);
        
        % transport fluxes    
        flux(Vex,ic) = TKinetics(model,pvec,mc(:,ic),Vex);
        
        % ETC fluxes
        flux = ETCflux(model,pvec,mc(:,ic),flux);

        % other fixed exchaged fluxes - currently sets them to 0  
        flux(VFex,ic) = EKinetics(model,pvec,mc(:,ic),VFex,flux(:,ic));
        
        % atp maintanance
        tatpm = strcmpi(model.rxns,'ATPM');
        sbid = model.S(:,tatpm)<0;
        sbid(h2o) = 0;
        flux(tatpm,ic) = pvec.Vmax(tatpm)*18.84*mc(sbid,ic)/pvec.K(sbid,tatpm)/...
                      (1+mc(sbid,ic)/pvec.K(sbid,tatpm));
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
        flux(model.bmrxn,ic) = model.Vss(model.bmrxn)/3600;
%         flux(strcmpi('GLCpts',model.rxns),ic) = 10;
    else
        %determine which group idx belongs to
        for id = 1:length(idx)
            if ismemeber(idx(id),Vind)
    %         if ismember(idx(id),vcup) || ismember(idx(id),vred)
    %             idflux(idx(id)) = flux(idx(id));
    %         elseif ismember(idx(id),Vind)
                idflux(idx(id),ic) = CKinetics(model,pvec,mc(:,ic),idx(id));
            elseif ismember(idx(id),Vex)
                idflux(idx(id),ic) = TKinetics(model,pvec,mc(:,ic),idx(id));
            elseif ismember(idx(id),VFex)
                idflux(idx(id),ic) = EKinetics(model,pvec,mc(:,ic),VFex,idflux);
            elseif idx(id)==find(strcmpi(model.rxns,'atpm'))
    %             if mc(strcmpi(model.mets,'atp[c]'))>0 &&...
    %                mc(strcmpi(model.mets,'h2o[c]'))>0
    %                 idflux(idx(id)) = 8.39;
    %             else
    %                 idflux(idx(id)) = 0;
    %             end
                idflux(idx(id),ic) = 0;
            elseif idx(id)==model.bmrxn
                idflux(id,ic) = 0;
            end        
        end        
        idxflux(strcmpi('GLCpts',model.rxns),ic) = 20;        
    end
end

if ~isempty(idx)
    flux = idxflux(idx,:);
end

% if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
%     all(mc(logical(model.S(:,model.bmrxn)>0))>0)
%     flux(model.bmrxn) = 0.1;
% else
%     flux(model.bmrxn) = 0;
% end

