%data-based elimimation of non-confirming models from an ensemble
%Obtain Metabolite data
    %Using Chassagnole's model instead
%Run network kinetic model for chosen enzyme and first perturbation
    [inSolution,enSolution] =  solveEnsemble(ensb,FBAmodel,variable,pertb,SolverOptions);
%Compare data points at different time points
%Inputs
nmodels = length(fieldnames(inSolution));
%1. Initial Perturbation
pertb.pertb1.enzname = 'Protein2';
pertb.pertb1.change = 'increase';
pertb.pertb1.percent = 0;
%2. Do perturbations
[enSolution,npertb] = doPerturbation(pertb,model,ensb,variable,inSolution,SolverOptions);
%3. Compare data points and eliminate models
%4. Iterate

%while condition is satisfied or not satisfied
%while accp_models 
    modelnames = fieldnames(ensb);
    nmodels = length(modelnames);
    %Choose perturbations
    % define #pertb
    % Solve Zommordi et al., 2013 optimization problem here?
    %Do Perturbation
    [enSolution,npertb] = doPerturbation(pertb,model,ensb,variable,inSolution,SolverOptions);
    %Compare data points
    %Choose Metabolite of i9nterest for comparison
    %define #metabolite name
    %define #time course concentrations (t,x)
    %define #time 
    %define #concentrations
    %[parsedmodels] = compareMetabolite(enSolution,nmodels);
    expt_t = data.t;
    expt_C = data.c;
    %tfm = strcmp(metabolite name, model.Metabolites);
    for imodel = 1:nmodels
        mname = modelnames{imodel};
        model_t = enSolution.pertb1.(mname).t;
        model_C = enSolution.pertb1.(mname).y;
        %Time comaprison
        for idata = 1:length(expt_t)
            sim_time = abs(model_t - expt_t(idata)) <= delt;
            if any(sim_time)
                simC = model_C(tfm,sim_time);
                avg_simC(idata) = sum(simC)/length(simC);               
            end
        end
         %Compare Data at time = expt_t(idata)
        if abs(avg_simC - expt_C) <= delC
            decision.(mname)(idata) = 1;
        end        
    end
    for imodel = 1:nmodels
        mname = modelnames{imodel};
        if all(decision.(mname))
            %keep model in ensemble
            new_ensb.(mname) = ensb.(mname);
            accp_models = accp_models + 1;
        else
            %Eliminate models that do not confirm to data points
            fprintf('Model %s has been eliminated\n',mname);
            elim_models = elim_models + 1;
        end
    end
    %Diagnostics
    elim_models;
    accp_models;
    %Reassignment
    ensb = new_ensb;
       
    
%end


