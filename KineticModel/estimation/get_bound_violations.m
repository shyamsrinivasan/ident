% compare optimization results with variable bounds
% similar to Laurence's work with Emilio
function get_bound_violations(optsol,prob,optimdata)

if isfield(prob,'lb')
    lb = prob.lb;
end
if isfield(prob,'ub')
    ub = prob.ub;
end
nvar = optimdata.nvar;
nc = optimdata.nc;
p_id = optimdata.p_id;
nf = optimdata.nf;
npert = optimdata.npert;

midpt = lb + (ub-lb)./2;
nc_midpt = midpt(1:nc*npert);
np_midpt = midpt(nc*npert+1:nc*npert+length(p_id));
nf_midpt = midpt(nvar-nf*npert+1:nvar);

% plot data
xdata_nc = 1:nc*npert;
xdata_np = 1:length(p_id);
xdata_nf = 1:nf*npert;

ydata_nc_1 = lb(1:nc*npert);
ydata_np_1 = lb(nc*npert+1:nc*npert+length(p_id));
ydata_nf_1 = lb(nvar-nf*npert+1:nvar);
ydata_nc_2 = ub(1:nc*npert);
ydata_np_2 = ub(nc*npert+1:nc*npert+length(p_id));
ydata_nf_2 = ub(nvar-nf*npert+1:nvar);

% think about replacing bar with line graphs with big line widths
% concentration figure
figure
bh = barh(xdata_nc,ydata_nc_1);
hold on % replace with axis property next figure
bh = barh(xdata_nc,ydata_nc_2);
% change base value with bar properties to midpt
% plot optsol values here

% parameter figure
figure
bh = barh(xdata_np,ydata_np_1);
hold on % replace with axis property next figure
bh = barh(xdata_np,ydata_np_2);
% change base value with bar properties
% plot optsol values here

% flux figure
figure
bh = barh(xdata_nf,ydata_nf_1);
hold on % replace with axis property next figure
bh = barh(xdata_nf,ydata_nf_2);
% change base value with bar properties
% plot optsol values here

