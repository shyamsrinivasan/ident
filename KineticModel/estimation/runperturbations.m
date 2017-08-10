function [xss_pval,fss_pval,collect_p,dyndata] = runperturbations(fh,pid,pval,opts)
if nargout>3
    getdyndata = 1;
end

np = length(pval);
xss_pval = zeros(length(opts.x0),np);
fss_pval = zeros(6,np);
odep_bkp = opts.odep;

n_odep = length(opts.odep);
collect_p = zeros(np,n_odep);

dyndata = struct();

for ival = 1:np
    opts.odep(pid) = pval(ival);
    [xdyn,fdyn,xss,fss] = fh(opts);
    xss_pval(:,ival) = xss;
    fss_pval(:,ival) = fss;
    if getdyndata
        dyndata(ival).xdyn = xdyn;
        dyndata(ival).fdyn = fdyn;
    end
    collect_p(ival,:) = opts.odep;
    opts.odep = odep_bkp;
end
