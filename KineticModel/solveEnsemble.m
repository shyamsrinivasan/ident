%function [inSolution,enSolution] =...
%          solveEnsemble(ensb,model,variable,pertb,SolverOptions)
%solve an ensemble of models 
%Need to store & plot flux data
function [inSolution,enSolution] =...
         solveEnsemble(ensb,model,variable,pertb)
nmodels = length(fieldnames(ensb));
saveData.dirname = 'C:\Users\shyam\Documents\Courses\CHE 1125 Project\Results\Kmodel';
%Initial Model Perturbation
fprintf('Initial Model Integration to obtain initial SS\n');
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    inSolution.(mname) = struct([]);
    fprintf('%s\n',mname);
%     saveData.filename = mname;
    %Model initialization
    [model,batch,solverP,saveData] = initializeModel(FBAmodel,10);
    %simulate models first to get initial SS
    Solution = callODEsolver(model,ensb.(mname),variable,inSolution.(mname),batch,solverP);
    inSolution.(mname) = Solution;  
    %Calculate initial & final flux   
    inSolution.(mname).flux = calc_flux(model,ensb.(mname),Solution.y(:,end));
    close all
end
%Enzyme Perturbations
fprintf('Model Perturbations\n');
%Function to perform different enzyme perturbations
[enSolution,npertb] = doPerturbation(pertb,model,ensb,variable,inSolution,SolverOptions);
% npertb = length(fieldnames(pertb));
% for isample = 1:nmodels
%     mname = sprintf('model%d',isample);
%     fprintf('%s\n',mname);
%     for ipertb = 1:npertb
%         pname = sprintf('pertb%d',ipertb);
%         fprintf('Perturbation %d\n',ipertb);
%         [Solution] = pertbEnzyme(pertb.(pname),model,ensb.(mname),...
%                                  variable,inSolution.(mname),SolverOptions);
%         enSolution.(pname).(mname) = Solution;        
%         enSolution.(pname).(mname).flux = calc_flux(model,ensb.(mname),Solution.y(:,end));
%     end
%     close all
% end
%save concentration and flux data for 1 model named "model1"
[nvar,ntpts] = size(enSolution.pertb1.model1.y);
save_data.t = enSolution.pertb1.model1.t;
save_flux.y = zeros(npertb,model.nt_rxn);
%save concentration and flux data for all models
for isample = 1:nmodels
    save_data.y = zeros(npertb*nvar,ntpts);
    mname = sprintf('model%d',isample);
    saveData.filename = mname;   
    save_flux.t = inSolution.(mname).flux;
    ipertb = 0;
    while ipertb < npertb
        pname = sprintf('pertb%d',ipertb+1);
        save_data.y(nvar*ipertb+1:nvar*(ipertb+1),:)=enSolution.(pname).(mname).y;
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
varname = {'A','B','C','D','E','P'};
[hfig,hsubfig] = compareEnsemble(model,npertb,varname,enSolution,ColorSpec,LineP);
%set default Line Properties
setProperties(hfig,hsubfig,enSolution.pertb1.model1);

%Plot - Comparison between fluxes
var_ind = [model.Vind,model.Vexind'];
F_ColorSpec = setLineP_flux(nmodels);
[h_ffig,h_fsubfig] = compareFluxes(nmodels,npertb,var_ind,enSolution,F_ColorSpec);
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


function setProperties(hfig,hsubfig,Solution)
if ~isempty(hsubfig)
    for ivar = 1:length(hsubfig)
        Xmax = max(Solution.t(:));
        Xdiv = sd_round(max(Solution.t(:))/10,1,1);
        %Common Axis Properties
        AxisP = struct();
%         AxisP.XMinorTick = 'on';
        AxisP.XTick = 0:Xdiv:Xmax;
%         AxisP.YMinorTick = 'on';
        AxisP.TickLength = [0.02,0.02];
        AxisP.XColor = [.4 .4 .4];
        AxisP.YColor = [.4 .4 .4];
        AxisP.XLim = [0 Xmax];
        AxisP.FontName = 'Courier';
        AxisP.FontSize = 14; 

        %Variable Specific Properties
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);
        hline = findobj(hca,'type','line');
        Y = zeros(length(hline),2);
        iline = 1;
        while iline <= length(hline)
            Y(iline,1) = min(get(hline(iline),'YData'));
            Y(iline,2) = max(get(hline(iline),'YData'));
            iline = iline+1;
        end  
%         Ymin = min(Y(:,1));
%         Ymax = max(Y(:,2));
%         Ymin = sd_round(min(Y(:,1)),3,5);
%         Ymax = sd_round(max(Y(:,2)),3,5);
%         Ydiv = max(Y(:,2))-min(Y(:,1))/10;
%         Ydiv = sd_round((max(Y(:,2))-min(Y(:,1)))/10,2,5);
%         AxisP.YLim = [Ymin Ymax];           
%         AxisP.YTick = Ymin:Ydiv:Ymax;
        if ~isempty(AxisP)
            set(hca,AxisP);
        end
    end
    legend(findobj(hsubfig(1),'type','axes'),'show','Location','Best');
    legend(findobj(hsubfig(1),'type','axes'),'boxoff');
end 
   