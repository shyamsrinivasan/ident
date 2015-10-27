function flux = iflux(model,pvec,mc,flux,idx)
if nargin <5
    idx = [];
else
    idflux = zeros(length(idx),1);
end
if nargin<4
    flux = zeros(model.nt_rxn,1);
end

%carbon uptake fluxes
[flux,~,vcup] = CarbonKinetics(model,pvec,mc,flux);

%redox reaction fluxes in vrem
[flux,~,vred] = RedoxKinetics(model,pvec,mc,flux);

%intracellular fluxes in Vind
Vind = model.Vind;
vrem = [find(strcmpi(model.rxns,'GLCpts'))...        
        find(strcmpi(model.rxns,'ATPM'))...
        vred];
Vind = setdiff(Vind,vrem);    

%transport fluxes
Vex = model.Vex;
Vex = setdiff(Vex,vrem);

%other fixed exchaged fluxes
VFex = model.VFex;

if isempty(idx)
    
    flux(Vind) = CKinetics(model,pvec,mc,Vind);

    %transport fluxes    
    flux(Vex) = TKinetics(model,pvec,mc,Vex);

    %other fixed exchaged fluxes    
    flux(VFex) = EKinetics(model,pvec,mc,VFex,flux);

    %biomass
    if mc(strcmpi(model.mets,'atp[c]'))>0 &&...
       mc(strcmpi(model.mets,'h2o[c]'))>0
        flux(strcmpi(model.rxns,'atpm')) = 8.39;
    else
        flux(strcmpi(model.rxns,'atpm')) = 0;
    end

    % flux(strcmpi(model.rxns,'atpm')) = 0;
    if mc(logical(model.S(:,model.bmrxn)<0))>0
        flux(model.bmrxn) = model.Vss(model.bmrxn);%0.01;
    else
        flux(model.bmrxn) = 0;
    end
else
    %determine which group idx belongs to
    for id = 1:length(idx)
        if ismember(idx(id),vcup) || ismember(idx(id),vred)
            idflux(idx(id)) = flux(idx(id));
        elseif ismember(idx(id),Vind)
            idflux(idx(id)) = CKinetics(model,pvec,mc,idx(id));
        elseif ismember(idx(id),Vex)
            idflux(idx(id)) = TKinetics(model,pvec,mc,idx(id));
        elseif ismember(idx(id),VFex)
            idflux(idx(id)) = EKinetics(model,pvec,mc,VFex,idflux);
        elseif idx(id)==find(strcmpi(model.rxns,'atpm'))
            if mc(strcmpi(model.mets,'atp[c]'))>0 &&...
               mc(strcmpi(model.mets,'h2o[c]'))>0
                idflux(idx(id)) = 8.39;
            else
                idflux(idx(id)) = 0;
            end
        elseif idx(id)==model.bmrxn
            idflux(id) = 0;
        end        
    end
    flux = idflux(idx);
end

% if flux(strcmpi(model.rxns,'atpm'))>=1e-5 &&...
%     all(mc(logical(model.S(:,model.bmrxn)>0))>0)
%     flux(model.bmrxn) = 0.1;
% else
%     flux(model.bmrxn) = 0;
% end

