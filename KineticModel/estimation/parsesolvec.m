function [xconc,xpar,xflx] = parsesolvec(xval,data)

np = length(data.p_id);
xconc = xval(1:data.nc*data.npert,:);
xpar = xval(data.nc*data.npert+1:data.nc*data.npert+np,:);
xflx = xval(data.nvar-data.nf*data.npert+1:data.nvar,:);

% % rearrange xconc to nc x npert matrix
% xconc = reshape(xconc,[data.nc data.npert]);
% % rearrange xflxeq to nf x npert matrix
% xflx = xflx';
