function x_opt = runoptimp(opts,plist,odep_bkp,optimstruct)

noptim = size(optimstruct,2);
x_opt = cell(noptim,1);

for j = 1:noptim
    if isfield(optimstruct,'xss')
        xss = optimstruct(j).xss;
        fss = optimstruct(j).fss;
    end
    x_opt{j} = optimize_p(opts,xss,fss,plist,odep_bkp);
end