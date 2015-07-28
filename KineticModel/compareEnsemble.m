%function [hfig,hsubfig] = compareEnsemble(model,npertb,varname,solution,...
%                                          ColorSpec,LineP,hsubfig)
%**************************************************************************
%Compare dynamic response of metabolites to input perturbations
%December 2014
%**************************************************************************
function [hfig,hsubfig] = compareEnsemble(model,npertb,varname,solution,...
                                          ColorSpec,LineP,hsubfig)
nvar = length(varname);
if nargin < 7
    hsubfig = zeros(nvar*npertb,1);
end
if nargin < 6
    LineP = struct();
end
y_label = cell(length(varname),1);
for ivar = 1:length(varname)
    var_tf = strcmpi(varname{ivar},model.mets);
    y_label{ivar} = sprintf('%s mmole/g DCW',model.mets{var_tf});
    if any(var_tf(1:model.nint_metab))            
        var_tf = logical(var_tf(1:model.nint_metab));
    end
    %Figure Properties, Misc
    if isempty(findobj('type','figure','Name','Compare Solutions'))
        hfig = figure('Name','Compare Solutions'); 
        figure(hfig);
    end     
    for ipertb = 1:npertb
        pname = sprintf('pertb%d',ipertb);
        iplot = ipertb+npertb*(ivar-1);
        if hsubfig(iplot) ~= 0
            hca = findobj(hsubfig(iplot),'type','axes');
            set(hfig,'CurrentAxes',hca);                
        else%subplot is unassigned
            hsubfig(iplot) = subplot(nvar,npertb,iplot);  
            hca = gca;
            %Make sure more plots can be added at end of loop
            set(hca,'NextPlot','add');     
            set(hsubfig(iplot),'NextPlot','add');
        end
        for imodel = 1:length(fieldnames(solution.(pname)))
            mname = sprintf('model%d',imodel);
            data = solution.(pname).(mname);
            LineP.DisplayName = mname;
            LineP.Color = ColorSpec{imodel};
            hline = line('XData',data.t(:),'YData',data.y(var_tf,:));            
            if ~isempty(LineP)
                set(hline,LineP);
            else
                set(hline,'Color',[0 .5 0]);
            end
        end
    end
end
%Identifying last row
last_hsub = reshape(hsubfig,[npertb,nvar])';
lr_hsub = last_hsub(end,:)';
for isub = 1:length(lr_hsub)
    x_label = 'Time(s)';   
    hca = findobj(lr_hsub(isub),'type','axes');
    set(hfig,'CurrentAxes',hca);
    set(get(gca,'XLabel'),'String',x_label);
    set(get(gca,'XLabel'),'FontName','Arabic Type Setting');
    set(get(gca,'XLabel'),'FontSize',12);
end
%Identifying first column
fc_hsub = last_hsub(:,1);
for isub = 1:length(fc_hsub)
    %x_label = 'Time(s)';   
    hca = findobj(fc_hsub(isub),'type','axes');
    set(hfig,'CurrentAxes',hca);
    set(get(gca,'YLabel'),'String',y_label{isub});  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'YLabel'),'FontSize',12);
end    
legend(findobj(fc_hsub(1),'type','axes'),'show');
legend(findobj(fc_hsub(1),'type','axes'),'boxoff');
%Identifying First row
fr_hsub = last_hsub(1,:);
for isub = 1:length(fr_hsub)
    hca = findobj(fr_hsub(isub),'type','axes');
    set(hfig,'CurrentAxes',hca);
    pname = sprintf('Perturbation %d',isub);
    set(gca,'Title',text('String',pname,'Color','k'));
    set(get(gca,'Title'),'FontName','Arabic Type Setting');
    set(get(gca,'Title'),'FontSize',12);
    set(get(gca,'Title'),'FontWeight','bold');
end
FProperty.NumberTitle = 'off';
FProperty.Color = [1 1 1];  
set(hfig,FProperty);
return
    