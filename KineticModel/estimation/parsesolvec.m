function [xconc,xpar] = parsesolvec(optsol,data)

np = length(data.p_id);
xconc = optsol.xval(1:data.nc);
xpar = optsol.xval(data.nc+1:data.nc+np);
