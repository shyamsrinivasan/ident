function sol = getperturbations(ptopts,fh,opts,sol,getdyndata)
if nargin<5
    getdyndata = 0;
end
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
    if ~getdyndata
        [xss,fss,collect_p] = runperturbations(fh,exp_pid,exp_pval,opts);
    else
        [xss,fss,collect_p,dyndata] = runperturbations(fh,exp_pid,exp_pval,opts);
    end
    sol(istart).xss = xss;
    sol(istart).fss = fss;
    if getdyndata
        sol(istart).xdyn = dyndata.xdyn;
        sol(istart).fdyn = dyndata.fdyn;
    end
    sol(istart).exp_pid = exp_pid;
    sol(istart).exp_pval = exp_pval;
    sol(istart).odep = collect_p;
    
    istart = istart+1;
end

