%Plot Both TRN and Kinetic Models (Time Course)
function [hfig,data,hsubfig] =...
         DynPlot(model,ng,varname,data,hfig,LineP,hsubfig,call)
if nargin < 8
    call = 0;
end
if isempty(ng)
    mdes = 2;%Kientic Model
else
    mdes = 1;%TRN Model    
end
if isfield(model,'nvar')
    nvar = model.nvar;
else
    nvar = size(data.y,1);
end
nadd = 0;
Var = cell(nvar,1);
if mdes == 1
    protind = ng(1)+(1:ng(2));
    cmind = ng(1)+model.PMind_R;
    mt_ind = ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)-ng(6)+1;
    mx_ind = sum(ng)-ng(6)+2:sum(ng)+1;
    Var(1:ng(1)) = model.Gene;
    Var(protind) = model.regs(protind-ng(1));
    Var(ng(1)+model.rnap_ind) = model.regs(model.rnap_ind);
    Var(cmind) = model.regs(cmind-ng(1));
    Var(mt_ind) = model.regs(mt_ind-ng(1));
    Var(mx_ind) = model.regs(mx_ind-ng(1));
    
    %data.t = data.t/3600;
    
    %Add protein name
    for ivar = 1:length(varname)
        var_tf = strcmpi(varname{ivar},Var);
        if any(var_tf(1:ng(1)))        
            nadd = nadd+1;
            varname{length(varname)+1} =...
            model.regs{logical(full(model.trate*var_tf(1:ng(1))))};        
        end    
    end
elseif mdes == 2
    Var = model.mets;
end

nplots = length(varname);
if nargin < 7 || isempty(hsubfig)
    hsubfig = zeros(nplots,1);
end

if nargin < 6 || isempty(LineP)
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
y_label1 = cell(nplots,1);
for ivar = 1:nplots
    var_tf = strcmpi(varname{ivar},Var);%Any Model
    if mdes == 1%TRN Model
        if any(var_tf(1:ng(1)))||any(var_tf(ng(1)+1:ng(1)+ng(2)))%||...
    %        any(var_tf(cx_ind))%mRNA or protein        
            g_tf = strcmpi(varname{ivar},model.Gene);
            pg_tf = strcmpi(varname{ivar},model.regs);
            if any(g_tf)
                %mRNA
                y_label1{ivar} = sprintf('%s mRNA umole/gDCW',model.Gene{g_tf});
                var_tf = logical([var_tf(1:ng(1));zeros(nvar-ng(1),1)]);     
    %             LineP.Displayname = sprintf('%s',model.Gene{g_tf});          
            elseif any(pg_tf)
                %Protein            
                y_label1{ivar} = sprintf('%s umole/gDCW',model.regs{pg_tf});
                var_tf = logical([zeros(ng(1),1);pg_tf;zeros(nvar-ng(1)-length(pg_tf),1)]);        
    %             LineP.Displayname = sprintf('%s',model.regs{pg_tf});            
            end
        elseif any(var_tf(ng(1)+ng(2)+1:ng(1)+ng(2)+ng(3)))
            cmp_tf = strcmpi(varname{ivar},model.regs);
            mc_tf = strcmpi(varname{ivar},model.mets);
            if any(cmp_tf)&&any(mc_tf)
                y_label1{ivar} = sprintf('%s umole/gDCW',model.regs{cmp_tf});
            end
        elseif any(var_tf(ng(1)+ng(2)+ng(3)+2:end))
            m_tf = strcmpi(varname{ivar},model.mets);        
    %         if any(m_tf) && var_tf(end)%biomass
    %             y_label1 = sprintf('%s (gDCW)',model.mets{end});
    %             LineP.Displayname = y_label1;
            if any(m_tf)%Metabolite
                y_label1{ivar} = sprintf('%s mmole/gDCW',model.mets{m_tf}); 
    %             LineP.Displayname = sprintf('%s',model.mets{m_tf});
            end      
        else
            fprintf('Variable %s does not Exist\n',varname{ivar});
            continue
        end
        x_label = 'Time (h)';
    elseif mdes == 2%Kinetic Model
        y_label1{ivar} = sprintf('%s mmole/g DCW',model.mets{var_tf});
        if any(var_tf(1:model.nint_metab))            
            var_tf = logical(var_tf(1:model.nint_metab));
        end
        x_label = 'Time (s)';
    end
    
    %Figure Properties, Data, etc.,
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
    hline1 = line('Xdata',data.t(:),'YData',data.y(var_tf,:));    
    if call %if called from trnperturbation
%         setProperties(hfig,hsubfig(ivar),hca,hline1,data.t(:),data.y(var_tf,:),LineP);
        setProperties(hfig,hsubfig,data);
        set(hline1,LineP);
    else
        set(hline1,LineP);
    end
    
    % whitebg(hfig,[0 0 0]);
    set(get(hca,'XLabel'),'String',x_label); 
    set(get(hca,'XLabel'),'FontName','Lucida Sans');
    set(get(hca,'XLabel'),'FontSize',10);
    set(get(hca,'YLabel'),'String',y_label1{ivar}); 
    set(get(hca,'YLabel'),'FontName','Lucida Sans');   
    set(get(hca,'YLabel'),'FontSize',10);
end

return