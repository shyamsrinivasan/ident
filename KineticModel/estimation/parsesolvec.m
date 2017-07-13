function [xconc,xpar] = parsesolvec(optsol,data)

np = data.nvar-data.nc;
xconc = optsol.xval(1:data.nc);
xpar = optsol.xval(data.nc+1:data.nc+np);

