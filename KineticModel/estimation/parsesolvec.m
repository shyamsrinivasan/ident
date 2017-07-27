function [xconc,xpar,xflx] = parsesolvec(optsol,data)

np = length(data.p_id);
xconc = optsol.xval(1:data.nc*data.npert);
xpar = optsol.xval(data.nc*data.npert+1:data.nc*data.npert+np);
xflx = optsol.xval(data.nvar-data.nf*data.npert+1:data.nvar);

% rearrange xconc to nc x npert matrix
xconc = reshape(xconc,[data.nc data.npert]);
% rearrange xflxeq to nf x npert matrix
xflx = xflx';
