% change lb and ub for specific parametyers for specific fluxes
function [new_lb,new_ub] = boundsf3(lb,ub,data)

nc = data.nc;
npert = data.npert;
np = data.np;

% parameter for flux 3
new_lb = lb;
new_ub = ub;

% set separate bounds for Kmfdp for flux 3
new_lb(nc*npert+1:nc*npert+1) = .03*ones(1,1); 
new_ub(nc*npert+1:nc*npert+1) = 1*ones(1,1);
% set separate bounds for Kmpep (reg) for flux 3
new_lb(nc*npert+2:nc*npert+2) = .03*ones(1,1); 
new_ub(nc*npert+2:nc*npert+2) = 1*ones(1,1);
% set separate bounds for kcat for flux 3
new_lb(nc*npert+3:nc*npert+np) = .1*ones(1,1); 
new_ub(nc*npert+3:nc*npert+np) = 2*ones(1,1);

