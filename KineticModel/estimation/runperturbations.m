function [xss_pval,fss_pval,collect_p] = runperturbations(fh,pid,pval,opts)

np = length(pval);
xss_pval = zeros(length(opts.x0),np);
fss_pval = zeros(6,np);
odep_bkp = opts.odep;

n_odep = length(opts.odep);
collect_p = zeros(np,n_odep);

for ival = 1:np
    opts.odep(pid) = pval(ival);
    [~,~,xss,fss] = fh(opts);
    xss_pval(:,ival) = xss;
    fss_pval(:,ival) = fss;
    collect_p(ival,:) = opts.odep;
    opts.odep = odep_bkp;
end