% parse given solution vector into concentrations, fluxes and parameters
function new_optsol = parse_optsol(optsol,data)

nval = size(optsol,2);
new_optsol = optsol;
if nval>1
    xval = cat(2,optsol.xval);
    [xconc,xpar,xflx,xnoise] = parsesolvec(xval,data);
else
    [xconc,xpar,xflx,xnoise] = parsesolvec(optsol.xval,data);
end

% reshape and add to structure array
for ival = 1:nval
    new_optsol(ival).xconc = reshape(xconc(:,ival),[data.nc,data.npert]);
    new_optsol(ival).xpar = xpar(:,ival);
    new_optsol(ival).xflx = xflx(:,ival)';
    if ~isempty(xnoise)
        new_optsol(ival).xnoise = xnoise(:,ival);
    end
end