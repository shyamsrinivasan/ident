function [h_fig,h_subfig] =...
         printMetResults(FBAmodel,allSolution,conc,petconc,flux,varname)
%Plot Solution Curves
LineP = struct();
LineP.LineWidth = 2.0;
nsampl = length(fieldnames(allSolution));

for ism = 1:nsampl
    if ism == 1
        sim_name = 'init';
        LineP.Color = rgb('Green');
        LineP.DisplayName = sprintf('Initial');
        LineP.LineStyle = '--';
    elseif ism > 1
        %Select model
        sim_name = sprintf('sample_%d',ism-1); 
        %Plot solution curve
        LineP.Color = 'Black';
        LineP.LineStyle = '-';
        LineP.Displayname = sim_name;
    end
    if exist('h_subfig','var')
        [h_fig,h_subfig] =...
        plotSims(sim_name,FBAmodel,varname,allSolution,LineP,h_subfig);
    else
        [h_fig,h_subfig] =...
        plotSims(sim_name,FBAmodel,varname,allSolution,LineP);
    end        
end
setProperties(h_fig,h_subfig,allSolution.(sim_name));
plotData = [];

return

function [h_fig,h_subfig] =...
         plotSims(sim_name,FBAmodel,varname,allSolution,LineP,h_subfig)
%Plot solution curve                                                
if isempty(findobj('type','figure','Name','Compare Solutions'))
    h_fig = figure('Name','Compare Solutions','Color',[1 1 1]);
else
    h_fig = findobj('type','figure','Name','Compare Solutions');
    figure(h_fig);
end 
if nargin > 6 || exist('h_subfig','var')
    [h_fig,~,h_subfig] =...
    DynPlot(FBAmodel,[],varname,allSolution.(sim_name),h_fig,LineP,h_subfig);
%     [h_fig,h_subfig] =...
%     plotMetabolites(FBAmodel,varname,allSolution.(sim_name),h_fig,LineP,h_subfig);
else
    [h_fig,~,h_subfig] =...
    DynPlot(FBAmodel,[],varname,allSolution.(sim_name),h_fig,LineP);
%     [h_fig,h_subfig] =...
%     plotMetabolites(FBAmodel,varname,allSolution.(sim_name),h_fig,LineP);
end

return
