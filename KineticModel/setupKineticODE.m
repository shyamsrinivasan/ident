function [newmodel,newpvec,Nimc,solverP,flux,dXdt,jacobian] =...
        setupKineticODE(model,pvec,mc,essrxn,Vupstruct,tmax)
% setup kinetic model for integration
if nargin<6
    tmax = 500;
end
if nargin<5
    Vupstruct = ([]);
end
if nargin<4
    essrxn = {};
end
% check if concentrations are initialized
if nargin < 3
    error('getiest:NoA','No initial concentrations');
end
if nargin<2
    error('getiest:NoA','No parameter vector');
else
    nmodels = size(pvec,2);
end

% initialize solver properties
[model,solverP] = imodel(model,'ode',tmax,essrxn,Vupstruct);

% remove water and protons (held constant) from consideration in the model
% integration phase
bkupmodel = model;
newpvec(nmodels) = struct();
for im = 1:nmodels
    [model,outpvec,newmc] = addremoveMets(bkupmodel,{'h2o[c]','h2o[e]'},pvec(im),mc);
    [newmodel,newoutpvec,nnewmc,cnstmet] =...
    remove_eMets(model,outpvec,newmc,[model.Vind model.Vex],...
    {'ac[e]','pyr[e]','g6p[e]','succ[e]','mal[e]','akg[e]','for[e]','fum[e]','h[e]','h[c]'});
    newpvec = copystruct(newoutpvec,newpvec,im);     
end
                       
% preinitialize Vind and Vex to reduce integration overhead in iflux.m                       
newmodel.Vind = addToVind(newmodel,newmodel.Vind,newmodel.rxn_add,newmodel.rxn_excep);
newmodel.rxn_excep = union(newmodel.rxn_excep,newmodel.rxns(newmodel.Vind));
newmodel.Vex = addToVind(newmodel,newmodel.Vex,[],newmodel.rxn_excep);  

% only initialize for varmets   
nvar = length(newmodel.mets)-length(find(cnstmet));
imc = zeros(nvar,1);
imc(imc==0) = 1;

% separate variables from constants and assign constants to a field in
% newmodel
Nimc = nnewmc(1:nvar);
PM = nnewmc(nvar+1:end);
Nimc(imc==0) = 0;
newmodel.imc = imc;
newmodel.imc(newmodel.imc==0) = 1;
newmodel.PM = PM;

% system check
flux = zeros(newmodel.nt_rxn,nmodels);
dXdt = zeros(nvar,nmodels);
jacobian = struct();
for im = 1:nmodels
    % calculate initial flux
    flux(:,im) = iflux(newmodel,newpvec(im),[Nimc.*imc;PM]);
    % toy model
    dXdt(:,im) = ToyODEmodel(0,Nimc,[],newmodel,newpvec(im));
    % get jacobian and eigen values and eigne vectors
    % [J,lambda,w] = getjacobian(Nimc,newpvec(im),newmodel);
%     jacobian(im).J = J;
%     jacobian(im).lambda = lambda;
%     jacobian(im).w = w;
end
end

function outstruct = copystruct(instruct,outstruct,id)
if isfield(instruct,'K')
    outstruct(id).K = instruct.K;     
end
if isfield(instruct,'Klb')
    outstruct(id).Klb = instruct.Klb;
end
if isfield(instruct,'Kub')
    outstruct(id).Kub = instruct.Kub;
end
if isfield(instruct,'KIact')
    outstruct(id).KIact = instruct.KIact;
end
if isfield(instruct,'KIihb')
    outstruct(id).KIihb = instruct.KIihb;
end
if isfield(instruct,'Vmax')
    outstruct(id).Vmax = instruct.Vmax;
end
if isfield(instruct,'kfwd')
    outstruct(id).kfwd = instruct.kfwd;
end
if isfield(instruct,'krev')
    outstruct(id).krev = instruct.krev;
end
if isfield(instruct,'delGr')
    outstruct(id).delGr = instruct.delGr;
end
if isfield(instruct,'Kin')
    outstruct(id).Kin = instruct.Kin;
end
if isfield(instruct,'Kind')
    outstruct(id).Kind = instruct.Kind;
end
if isfield(instruct,'feasible')
    outstruct(id).feasible = instruct.feasible;
end     
end
