function [xconc,xpar,xflx,xnoise] = parsesolvec(xval,data)

np = data.np;
npert = data.npert;
nvar = data.nvar;
nc = data.nc;
nf = data.nf;
if isfield(data,'type')
    type = data.type;
else
    type = 1;
end

xconc = xval(1:nc*npert,:);
xpar = xval(nc*npert+1:nc*npert+np,:);
xflx = xval(nc*npert+np+1:nc*npert+np+nf*npert,:);
if type==2
    xnoise = xval(nvar-2:nvar);
else
    xnoise = [];
end

% % rearrange xconc to nc x npert matrix
% xconc = reshape(xconc,[data.nc data.npert]);
% % rearrange xflxeq to nf x npert matrix
% xflx = xflx';
