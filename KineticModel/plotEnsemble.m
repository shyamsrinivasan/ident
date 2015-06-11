%function [h_vec] = dynamicplot(trnmodel,gname,pname,data,figh)
%**************************************************************************
%Plot dynamic response of genes and proteins to input perturbations
%December 2013
%**************************************************************************
function [hsubfig] = plotEnsemble(model,ng,varname,data,hfig,LineP,hsubfig)
if nargin < 9
    savefile = struct([]);
end
if ~isempty(varname)
    nplots = length(varname);
end
if nargin < 7
    hsubfig = zeros(nplots,1);
end
if nargin < 6 
    LineP = struct();
    LineP.LineWidth = 1.5;
    LineP.Color = [0 .5 0];
end
if nargin < 5
    hfig = figure;    
else
    figure(hfig);    
end
if rem(nplots,2) == 0  
   n = nplots/2;
else
    n = (nplots+1)/2;
end
for ivar = 1:nplots
    %Data for plotting 
    if ng(1) == 0%Metabolic Model
        var_tf = strcmpi(varname{ivar},model.Metabolites);
        y_label = sprintf('%s Conentration',model.Metabolites{var_tf});
        if any(var_tf(1:model.nint_metab))            
            var_tf = logical(var_tf(1:model.nint_metab));
        end
    end
    %Figure Properties, Misc
    figure(hfig);
    if hsubfig(ivar) ~= 0
        hca = findobj(hsubfig(ivar),'type','axes');
        set(hfig,'CurrentAxes',hca);                
    else%subplot is unassigned
        hsubfig(ivar) = subplot(2,n,ivar);       
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
    x_label = 'Time(s)';    
    set(get(gca,'XLabel'),'String',x_label);
    set(get(gca,'XLabel'),'FontName','Arabic Type Setting');
    set(get(gca,'XLabel'),'FontSize',12);
    set(get(gca,'YLabel'),'String',y_label);  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'YLabel'),'FontSize',12);
    % whitebg(hfig,[0 0 0]);    
end        

end
