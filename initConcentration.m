
% [Y,Yflux] =...
% initializeConcentration(model,pmeter,variable,Vname,Vconc,type,initSol)
% Formulate the initial conecntration and flux vector from given initial 
% concentrations and input concentrations
function [Y,batch] =...
initConcentration(model,variable,type,Vname,Vconc,initSol)
if nargin < 6
    initSol = {};
end
if nargin < 5
    Vconc = [2000;1e-4;10000;1e-4;...
             0;100;1000];%mmoles 
end
if nargin < 4
    Vname = {'glc[e]';'h[e]';'h2o[e]';'pi[e]';...             
             'nh4[e]';'co2[e]';'o2[e]'};
end
if nargin < 3
    type = 1;
end
if nargin < 2
    variable.MC = zeros(model.nt_metab,1);
end
if type == 1%Kinetic Model
    nint_metab = model.nint_metab;
    next_metab = model.next_metab;
    nvar = model.nt_metab;%All Metabolites
    Y = zeros(nvar,1);
    if ~isemptyr(initSol)
        Y = initSol.y(:,end);
%         Yflux = initSol.flux;
    else
%         Mbio = strcmpi('biomass',model.mets);

%         Y(1:(nint_metab-1)) = variable.MC(1:(nint_metab-1));
        Y(1:nint_metab) = variable.MC(1:nint_metab);

%         Y(Mbio) = 0;
        Y(nint_metab+1:nint_metab+next_metab) = ...
        assign_extconc(Vname,Vconc,model);
%         Yflux = calc_flux(model,pmeter,Y);
%         Yflux = [];
    end
elseif type == 2%TRN Model
end
batch.init{1} = Vname;
batch.init{2} = Vconc;