function [xss_pval,fss_pval] = runperturbations(fh,pid,pval,opts)

np = length(pval);
xss_pval = zeros(length(opts.x0),np);
fss_pval = zeros(5,np);
% odep_bkp = opts.odep;

for ival = 1:np
    opts.odep(pid) = pval(ival);
    [~,~,xss,fss] = fh(opts);
    xss_pval(:,ival) = xss;
    fss_pval(:,ival) = fss;
end