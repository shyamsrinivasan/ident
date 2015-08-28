% Formulate the initial flux vector from given initial 
% concentrations and input concentrations
function Yflux =...
initFlux(model,pmeter,Y,type,initSol)
if nargin < 5
    initSol = {};
end
if type == 1%Kinetic Model
%     nint_metab = model.nint_metab;
%     next_metab = model.next_metab;
%     nvar = model.nt_metab;%All Metabolites
%     Y = zeros(nvar,1);
    if ~isemptyr(initSol)
%         Y = initSol.y(:,end);
        Yflux = initSol.flux;
    else
%         Mbio = strcmpi('biomass',model.mets);

%         Y(1:(nint_metab-1)) = variable.MC(1:(nint_metab-1));
%         Y(1:nint_metab) = variable.MC(1:nint_metab);

%         Y(Mbio) = 0;
%         Y(nint_metab+1:nint_metab+next_metab) = ...
%         assign_extconc(Vname,Vconc,model);
        Yflux = calc_flux(model,pmeter,Y);        
    end
elseif type == 2%TRN Model
end