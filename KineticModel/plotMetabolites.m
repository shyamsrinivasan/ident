function [hfig,hsubfig] =...
         plotMetabolites(model,varname,data,hfig,LineP,hsubfig)

nvar = length(varname);
if nargin < 6
    hsubfig = zeros(nvar,1);
end
if nargin < 5
    LineP = struct();
end
y_label = cell(length(varname),1);
for ivar = 1:length(varname)
    var_tf = strcmpi(varname{ivar},model.Metabolites);
    y_label{ivar} = sprintf('%s mmole/g DCW',model.Metabolites{var_tf});
    if any(var_tf(1:model.nint_metab))            
        var_tf = logical(var_tf(1:model.nint_metab));
    end
    %Figure Properties, Misc
    if isempty(findobj('type','figure','Name','Compare Solutions'))
        hfig = figure('Name','Compare Solutions'); 
        figure(hfig);
    end   
    if hsubfig(ivar) ~= 0
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(ivar) = subplot(2,nvar,ivar);  
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(ivar),'NextPlot','add');
    end   
    hline = line('XData',data.t(:),'YData',data.y(var_tf,:));
    if ~isempty(LineP)
        set(hline,LineP);
    else
        set(hline,'Color',[0 .5 0]);
    end
end
return