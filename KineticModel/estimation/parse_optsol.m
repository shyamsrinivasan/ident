% parse given solution vector into concentrations, fluxes and parameters
function new_optsol = parse_optsol(optsol,data)

nval = size(optsol,2);
new_optsol = optsol;
if isfield(data,'pscale')
    pscale = data.pscale;
else
    pscale = [];
end
if nval>1
    xval = cat(2,optsol.xval);
    [xconc,xpar,xflx,xnoise] = parsesolvec(xval,data);
else
    [xconc,xpar,xflx,xnoise] = parsesolvec(optsol.xval,data);
end

% reshape and add to structure array
for ival = 1:nval
    new_optsol(ival).xss = reshape(xconc(:,ival),[data.nc,data.npert]);
    if ~isempty(pscale)
        new_optsol(ival).xpar = xpar(:,ival)./pscale;
    else
        new_optsol(ival).xpar = xpar(:,ival);
    end
    new_optsol(ival).fss = reshape(xflx(:,ival),[data.nf,data.npert]);
    if ~isempty(xnoise)
        new_optsol(ival).eps_cf = xnoise(:,ival);
    end
end