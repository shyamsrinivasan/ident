function varargout =...
solveEnsembleMC(model,ensb,variable,inSolution,varname,type)

if nargin < 6
    type = '';
end
nmodels = length(fieldnames(ensb));
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kmodel';

%Initial Model Perturbation
fprintf('Initial Model Integration to obtain initial SS\n');
for imodel = 1:nmodels
    mname = sprintf('model%d',imodel);
    if isempty(inSolution)
        inSolution.(mname) = struct([]);
    end
    fprintf('%s\n',mname);
%     saveData.filename = mname;
    %Model initialization
    [model,batch,solverP,saveData] = initializeModel(model,100);
    if nmodels > 1
        if isempty(inSolution.(mname))
            %simulate models first to get initial SS
            [Solution,finalSS] =...
            callODEsolver(model,ensb.(mname),variable,inSolution.(mname),batch,solverP);
            inSolution.(mname) = Solution;  
            
            %Calculate initial & final flux   
            inSolution.(mname).flux = calc_flux(model,ensb.(mname),Solution.y(:,end));
            allSolution.(mname).init = Solution;
            allfinalSS.(mname).init = finalSS;
            allfinalSS.(mname).init.flux = inSolution.(mname).flux;
            
            %Plot Initial Solution
            [hfig,hsubfig] =...
            printMetResults(model,allSolution,[],[],[],varname);
        end
    else
        if isempty(inSolution.model1)
            %simulate models first to get initial SS
            [Solution,finalSS] =...
            callODEsolver(model,ensb.(mname),variable,inSolution,batch,solverP);
            inSolution = Solution;             
            
            %Calculate initial & final flux   
            inSolution.flux = calc_flux(model,ensb.(mname),Solution.y(:,end));
            allSolution.init = Solution;
            allfinalSS.init = finalSS;   
            allfinalSS.init.flux = inSolution.flux;
            
            %Plot Initial Solution
            [hfig,hsubfig] =...
            printMetResults(model,allSolution,[],[],[],varname);
        end
    end
    close all
end

%Call MCsmulation for inital value MC on FBAmodel
nsamples = 500; %# samples
lb = [1e-4];
ub = [1000];
pvar = {'C[c]'};

if strcmpi(type,'MC')
    for imodel = 1:nmodels
        mname = sprintf('model%d',imodel);
        if nmodels > 1
            [MCdyn,MCss] =...
            MCsimulation(model,ensb.(mname),variable,nsamples,pvar,lb,ub,...
                         allSolution.(mname),allfinalSS.(mname),varname);
            allfinalSS.(mname) = MCss;
            allSolution.(mname) = MCdyn;
        else
            if nsamples > 1
                %Parallel
                [MCdyn,MCss,y0new] =...
                 MCsimulation_parallel(model,ensb.model1,variable,nsamples,pvar,lb,ub,...
                                       allSolution,allfinalSS,varname);
            else
        %        Serial
                [MCdyn,MCss] =...
                MCsimulation(model,ensb.model1,variable,nsamples,pvar,lb,ub,...
                             allSolution,allfinalSS,varname);
            end
            allSolution = MCdyn;
            allfinalSS = MCss;
        end
    end
    varargout{1} = allSolution;
    varargout{2} = allfinalSS;
    varargout{3} = y0new;
else
%Call enzyme pertubation for enzyme perturbation on ensemble
%     if nmodels > 1
    fprintf('Model Perturbations\n');
    %Function to perform different enzyme perturbations
    [enSolution,npertb] =...
    doPerturbation(pertb,model,ensb,variable,inSolution,SolverOptions);

    %save concentration and flux data for 1 model named "model1"
    [nvar,ntpts] = size(enSolution.pertb1.model1.y);
    save_data.t = enSolution.pertb1.model1.t;
    save_flux.y = zeros(npertb,model.nt_rxn);
    %save concentration and flux data for all models
    for imodel = 1:nmodels
        save_data.y = zeros(npertb*nvar,ntpts);
        mname = sprintf('model%d',imodel);
        saveData.filename = mname;   
        save_flux.t = inSolution.(mname).flux;
        ipertb = 0;
        while ipertb < npertb
            pname = sprintf('pertb%d',ipertb+1);
            save_data.y(nvar*ipertb+1:nvar*(ipertb+1),:)=...
            enSolution.(pname).(mname).y;
            save_flux.y(ipertb+1,:) = enSolution.(pname).(mname).flux;
            ipertb = ipertb+1;
        end
        save_data.y = [inSolution.(mname).y;save_data.y];
        savefile(save_data,mname,saveData);
        fname = sprintf('flux_%s',mname);
        savefile(save_flux,fname,saveData);
        %write concentrations to excel file
        status = selectData_kmodel(save_data,model,mname);
        status = selectFlux_kmodel(save_flux,model,fname);
    end
    
    %Plots - Comparison between models
    LineP = struct();
    LineP.LineWidth = 1.5;
    ColorSpec = setLineP(nmodels);
    [hfig,hsubfig] =...
    compareEnsemble(model,npertb,varname,enSolution,ColorSpec,LineP);
    %set default Line Properties
    setProperties(hfig,hsubfig,enSolution.pertb1.model1);

    %Plot - Comparison between fluxes
    var_ind = [model.Vind,model.Vexind'];
    F_ColorSpec = setLineP_flux(nmodels);
    [h_ffig,h_fsubfig] =...
    compareFluxes(nmodels,npertb,var_ind,enSolution,F_ColorSpec);
    varargout{3} = enSolution;
%     end
end
return

function ColorSpec = setLineP_flux(nlines)
ColorSpec = cell(nlines,1);
for icolor = 1:nlines
    fac = icolor/(1+icolor);
    ColorSpec{icolor} = [1/2*fac 1/2*fac 1/2*fac];  
end

function ColorSpec = setLineP(nlines)
ColorSpec = cell(nlines,1);
for icolor = 1:nlines
    fac = icolor/(1+icolor);
    ColorSpec{icolor} = [1/2*fac 3/4*fac 2/3*fac];  
end


