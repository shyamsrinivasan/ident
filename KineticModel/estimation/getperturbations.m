function sol = getperturbations(ptopts,fh,opts,sol)
npt = size(ptopts,2);
if nargin<4
    sol = struct([]);
    istart = 1;
%     iend = npt;
else
    istart = size(sol,2)+1;
%     iend = size(sol,2)+npt;
end

for i = 1:npt
    if isfield(ptopts,'exp_pid')
        exp_pid = ptopts(i).exp_pid;
    end
    if isfield(ptopts,'exp_pval')
        exp_pval = ptopts(i).exp_pval;
    end
    [xss,fss,collect_p] = runperturbations(fh,exp_pid,exp_pval,opts);
    sol(istart).xss = xss;
    sol(istart).fss = fss;
    sol(istart).exp_pid = exp_pid;
    sol(istart).exp_pval = exp_pval;
    sol(istart).odep = collect_p;
    
    istart = istart+1;
end

