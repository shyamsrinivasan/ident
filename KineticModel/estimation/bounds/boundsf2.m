% change lb and ub for specific parametyers for specific fluxes
function [new_lb,new_ub] = boundsf2(lb,ub,data)

nc = data.nc;
npert = data.npert;
p_id = data.p_id;

% parameter for flux 2
new_lb = lb;
new_ub = ub;

% set separate bounds for Km for flux 1
new_lb(nc*npert+1:nc*npert+1) = .03*ones(1,1); 
new_ub(nc*npert+1:nc*npert+1) = 1*ones(1,1);
% set separate bounds for Vmax for flux 1
new_lb(nc*npert+2:nc*npert+length(p_id)) = .1*ones(1,1); 
new_ub(nc*npert+2:nc*npert+length(p_id)) = 2*ones(1,1);

