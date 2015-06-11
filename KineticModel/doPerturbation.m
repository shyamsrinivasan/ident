function [ensembleSol,npertb] = doPerturbation(pertb,model,ensb,variable,initSol,SolverOptions)

nmodels = length(fieldnames(ensb));
npertb = length(fieldnames(pertb));
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    fprintf('%s\n',mname);
    for ipertb = 1:npertb
        pname = sprintf('pertb%d',ipertb);
        fprintf('Perturbation %d\n',ipertb);
        [Solution] = pertbEnzyme(pertb.(pname),model,ensb.(mname),...
                                 variable,initSol.(mname),SolverOptions);
        ensembleSol.(pname).(mname) = Solution;        
        ensembleSol.(pname).(mname).flux = calc_flux(model,ensb.(mname),Solution.y(:,end));
    end
    close all
end
return