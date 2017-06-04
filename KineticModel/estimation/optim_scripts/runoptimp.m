function x_opt = runoptimp(opts,plist,odep_bkp,optimstruct,optimfh)

noptim = size(optimstruct,2);
x_opt = cell(noptim,1);

for j = 1:noptim
    if isfield(optimstruct,'xss')
        xss = optimstruct(j).xss;
        fss = optimstruct(j).fss;
    end
    x_opt{j} = optimfh(opts,xss,fss,plist,odep_bkp);
end