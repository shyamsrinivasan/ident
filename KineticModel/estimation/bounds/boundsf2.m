% change lb and ub for specific parametyers for specific fluxes
function [new_lb,new_ub] = boundsf2(lb,ub,data)

nc = data.nc;
npert = data.npert;
np = data.np;

% parameter for flux 2
new_lb = lb;
new_ub = ub;

% set separate bounds for Km for flux 1
new_lb(nc*npert+1:nc*npert+1) = .008*ones(1,1); 
new_ub(nc*npert+1:nc*npert+1) = 1*ones(1,1);
% set separate bounds for Vmax for flux 1
new_lb(nc*npert+2:nc*npert+np) = .05*ones(1,1); 
new_ub(nc*npert+2:nc*npert+np) = 3*ones(1,1);

