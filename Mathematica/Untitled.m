%Split fluxes into forward and reverse fluxes
%Generate upper and lower bounds for fluxes, delGr and x(concentrations)
%Form A, bl and bu

%From model
S=TMFAmodel.S;
vl=TMFAmodel.vl; vu=TMFAmodel.vu; xl=TMFAmodel.xl; xu=TMFAmodel.xu;
dGl=TMFAmodel.dGl; dGu=TMFAmodel.dGu; dGo=TMFAmodel.dGo;

%   TMFAmodel: (structure) containing...
%       rxns [n x 1] reaction abbreviations
%       mets [m x 1] metabolite abbreviations
%       c [1 x n] objective function on fluxes
%       cConc [1 x m] objective function on concentrations
%       cGibbs [1 x n] objective function on Gibbs free energies
%       b [m x 1] mass balance constraint vector
%       S [m x n] stoichiometric matrix
%       vl [n x 1] flux lower bounds
%       vu [n x 1] flux upper bounds
%       xl [m x 1] concentration lower bounds (mM)
%       xu [m x 1] concentration upper bounds (mM)
%       dGo [nKnown x 1] standard reaction Gibbs change (if known)
%       dGoknown[nKnown x 1] indices of reactions with known dGo
%       dGl [numKnown x 1] reaction Gibbs change lower bound
%       dGu [numKnown x 1] reaction Gibbs change upper bound
%       uncertainty [numKnown x 1] uncertainty in dGo
%       subSystems {n x 1} (cell) reaction's subsystem classification
%       T [scalar] temperature in Kelvin

n=size(S,2); m=size(S,1);%m -metabolites, n- rxns
vflength=n; vrlength=n;
lnxlength=m;
dGflength=n;
zflength=n; zrlength=n;
ulength=vflength+vrlength+lnxlength+dGflength+zflength+zrlength;

vfstart=0;
vrstart=vfstart+vflength;
lnxstart=vrstart+vrlength;
dGfstart=lnxstart+lnxlength;
zfstart=dGfstart+dGflength;
zrstart=zfstart+zflength;

pvf=sparse(1:vflength,vfstart+(1:vflength),1,vflength,ulength);
pvr=sparse(1:vrlength,vrstart+(1:vrlength),1,vrlength,ulength);
plnx=sparse(1:lnxlength,lnxstart+(1:lnxlength),1,lnxlength,ulength);
pdGf=sparse(1:dGflength,dGfstart+(1:dGflength),1,dGflength,ulength);
pzf=sparse(1:zflength,zfstart+(1:zflength),1,zflength,ulength);
pzr=sparse(1:zrlength,zrstart+(1:zrlength),1,zrlength,ulength);

pv=sparse([1:n 1:n],[1:n vrstart+(1:n)],[ones(1,n) -ones(1,n)],n,ulength);








