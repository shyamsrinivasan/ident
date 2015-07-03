%function [] = plotSSexpression(trnmodel,ssdata,gname,pname,ngene,...
%                               tnreg,saveData,allSolution)
%Notes:
%Plot dynamic response of genes and proteins to input perturbations
%December 2013
%Very similar to dynamicplot.m - Was created as an alternative
%January 2014 - Obselete/Repurposed
%Written as a function plotSSexpression.m
%January 8th 2014
%Serves to plot just SS expressionvalues for comparison purposes
%Plots both dynamic and SS expression values for comparison purposes
%Plot each case (initialSS, perturbation 1, etc) separately
%January 11th 2014
%the steady state comparison plots require SS data
%from each experiment to be plotted as separate figures. This will be fixed
%in a later iteration to mimic the case for dynamic information.
%January 14th 2014
%Need to support the ability to plot specific experiment(s) only.
%**************************************************************************
function plotSSexpression(model,ng,conc,petconc,cncind,variable,hsubfig)
if isempty(ng)
    mdes = 2;%Kinetic Model
else
    mdes = 1;%TRN Model
end
if nargin < 3 || (length(cncind)==1 && strcmpi('all',cncind{1}))
    cncind = model.Regulators;
    ncnc = length(cncind);
else
    ncnc = length(cncind);
end
if rem(ncnc,2)==0
    n = ncnc/2;
else
    n = (ncnc+1)/2;
end
nvar = size(conc,1);
if nargin < 7
    hsubfig = zeros(nvar,1);
end
%Bins for all concentyrations
%ssconc       %row - concentrations %column - samples
conct = conc';%row - samples %column - concentrations
petconct = petconc';
% conciqr = iqr(conct);
% maxconc = max(conct);
% minconc = min(conct);
% h = 2.*conciqr.*(nsmpl^(-1/3));
% nbins = (maxconc-minconc)./h;

%Randomly choose colors
ColorSpec = chooseColors(ncnc);

Var = cell(nvar,1);
switch variable
    case 'concentration'
        if mdes == 1
            protind = ng(1)+(1:ng(2));
            cmind = ng(1)+model.PMind_R;
            mt_ind = ng(1)+ng(2)+ng(3)+2:sum(ng)-ng(5)+1;
            mx_ind = sum(ng)-ng(5)+2:sum(ng)+1;

            Var(1:ng(1)) = model.Gene;
            % Var(ng(1)+1:end) = model.Regulators;
            Var(protind) = model.Regulators(protind-ng(1));
            Var(ng(1)+model.rnap_ind) = model.Regulators(model.rnap_ind);
            Var(cmind) = model.Regulators(cmind-ng(1));
            Var(mt_ind) = model.Regulators(mt_ind-ng(1));
            Var(mx_ind) = model.Regulators(mx_ind-ng(1));
        elseif mdes == 2
            Var = model.mets;
        end        
    case 'flux'
        if mdes == 1 && isfield(model,'Regulators')            
            Var = model.Regulators;
        elseif mdes == 2 && isfield(model,'rxns')
            Var = model.rxns;
        end
end

for j = 1:ncnc
    
    switch variable   
        %Concentrations
        case 'concentration'
            tf = strcmpi(cncind{j},Var);            
            if mdes == 1
                if any(tf(1:ng(1)))||any(tf(ng(1)+1:ng(1)+ng(2)))
                    g_tf = strcmpi(cncind{j},Var(1:ng(1)));
                    p_tf = strcmpi(cncind{j},Var(ng(1)+1:ng(1)+ng(2)));
                    %Gene
                    if any(g_tf)
                        y_label1 = sprintf('%s mRNA umole/gDCW',model.Gene{g_tf});
                        tf = logical([tf(1:ng(1));zeros(nvar-ng(1),1)]); 
                    %Protein    
                    elseif any(p_tf)
                        y_label1 = sprintf('%s umole/gDCW',model.Regulators{p_tf});
                        tf = logical([zeros(ng(1),1);...
                                      p_tf;...
                                      zeros(nvar-ng(1)-length(p_tf),1)]); 
                    end
                elseif any(tf(ng(1)+ng(2)+2:end))
                    m_tf = strcmpi(cncind{j},Var(ng(1)+ng(2)+2:end));
                    %Metabolite
                    if any(m_tf)
                        y_label1 = sprintf('%s mmole/gDCW',model.mets{m_tf}); 
                    end
                else
                    fprintf('Variable %s does not Exist\n',cncind{j});
                    continue
                end
            elseif mdes == 2
                tf = strcmpi(cncind{j},Var);
                %Metabolite
                if any(tf)
                    y_label1 = sprintf('%s mmole/gDCW',model.mets{tf});
                else
                    fprintf('Variable %s does not Exist\n',cncind{j});
                    continue
                end                
            end
            if isempty(findobj('type','figure','Name','Steady State Concentrations'))
                hfig = figure('Name','Steady State Concentrations'); 
    %             figure(hfig);
            else
                hfig = findobj('type','figure','Name','Steady State Concentrations');
    %             figure(hfig);
            end  
            figure(hfig);
            hsubfig = lineplot(j,n,hfig,hsubfig,conct,petconct,tf,y_label1,ColorSpec);
        %Fluxes
        case 'flux'
            tfp = find(strcmpi(cncind{j},Var));
            %Protein/Flux
            if any(tfp)
                y_label1 = sprintf('Flux %s',Var{tfp});
            else
                fprintf('Variable %s does not Exist\n',cncind{j});
                continue
            end
            if isempty(findobj('type','figure','Name','Steady State Fluxes'))
                hfig = figure('Name','Steady State Fluxes'); 
    %             figure(hfig);
            else
                hfig = findobj('type','figure','Name','Steady State Fluxes');
    %             figure(hfig);
            end
            figure(hfig);
            hsubfig = lineplot(j,n,hfig,hsubfig,conct,petconct,tfp,y_label1,ColorSpec);     
    end
end
        
%Old parts
% hsubfig = zeros(ncnc,1);
% for j = 1:ncnc
%     tf = strcmpi(cncind{j},Var);
%     if mdes == 1
%         if any(tf(1:ng(1)))||any(tf(ng(1)+1:ng(1)+ng(2)))
%             g_tf = strcmpi(cncind{j},model.Gene);
%             p_tf = strcmpi(cncind{j},model.Regulators);
%             if any(g_tf)
%                 %mRNA
%                 y_label1 = sprintf('%s mRNA umole/gDCW',model.Gene{g_tf});
%                 tf = logical([tf(1:ng(1));zeros(nvar-ng(1),1)]);     
%     %             LineP.Displayname = sprintf('%s',model.Gene{g_tf});          
%             elseif any(p_tf)
%                 %Protein            
%                 y_label1 = sprintf('%s umole/gDCW',model.Regulators{p_tf});
%                 tf = logical([zeros(ng(1),1);p_tf;zeros(nvar-ng(1)-length(p_tf),1)]);        
%     %             LineP.Displayname = sprintf('%s',model.Regulators{pg_tf});            
%             end
%         elseif any(tf(ng(1)+ng(2)+2:end))
%             m_tf = strcmpi(cncind{j},model.mets);
%             if any(m_tf)%Metabolite
%                 y_label1 = sprintf('%s mmole/gDCW',model.mets{m_tf}); 
%         %             LineP.Displayname = sprintf('%s',model.mets{m_tf});
%             end
%         else
%             fprintf('Variable %s does not Exist\n',cncind{j});
%             continue
%         end
%     elseif mdes == 2
%         mtf = strcmpi(cncind{j},model.mets);
%         if any(mtf)
%             y_label1 = sprintf('%s mmole/gDCW',model.mets{mtf});
%         end
%     end
%     if isempty(findobj('type','figure','Name','Steady State Concentrations'))
%         hfig = figure('Name','Steady State Concentrations'); 
%         figure(hfig);
%     else
%         hfig = findobj('type','figure','Name','Steady State Concentrations');
%         figure(hfig);
%     end 
%     figure(hfig);
%     hsubfig(j) = subplot(2,n,j);
%     %Plot only initial and final(sample) concentrations
%     concd = conct(:,tf);
%     pconcd = petconct(:,tf);
%     X = zeros(length(concd)+length(pconcd),1);        
%     initc = concd(1);
%     concd(1) = [];
%     X(1) = 1;
%     X(2:length(pconcd)+1) = 1;
%     X(length(pconcd)+2:length(concd)+length(pconcd)+1) = 2;
%     hline = line(X,[initc;pconcd;concd]);
%     set(hline,'LineStyle','-',...
%               'Marker','o',...
%               'MarkerSize',5,...
%               'MarkerFaceColor',ColorSpec{j},...
%               'MarkerEdgeColor',ColorSpec{j});
% %     scatter(X,[initc;pconcd;concd],'MarkerFillColor',ColorSpec{j},...
% %                                    'MarkerEdgeColor',ColorSpec{j}); 
%     set(get(gca,'YLabel'),'String',y_label1);  
%     set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
%     set(get(gca,'YLabel'),'FontSize',12);  
%     xlabel = sprintf('Steady States');
%     set(get(gca,'XLabel'),'String',xlabel);  
%     set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
%     set(get(gca,'XLabel'),'FontSize',12);
%     Ydata = get(hline,'YData');
%     maxY = max(Ydata);
%     minY = min(Ydata);
%     if minY >= 1e-4 && maxY >= 1e-4
%         if maxY-minY <= 1e-2
%             maxY = maxY+1e-1;
%             if minY-1e-1 <= 0
%                 minY = 0;
%             else
%                 minY = minY-1e-1;
%             end
%         end
%     else
%         if maxY-minY <= 1e-4
%         maxY = maxY+1e-4;
%             if minY-1e-4 <= 0
%                 minY = 0;
%             else
%                 minY = minY-1e-4;
%             end
%         end
%     end
%     set(gca,'YLim',[minY maxY]);
%     AP.XColor = [.1 .1 .1];
%     AP.YColor = [.1 .1 .1];
%     AP.XLim = [0 max(X)+1];
%     AP.XTick = 0:1:max(X);
%     set(hsubfig(j),AP);    
%    
%     %Old Parts end
%     %determine number of elements and indices of elements in each bin
% %     binEdge = linspace(max(concd),min(concd),nbins+1);
% %     binUEdge = binEdge(2:end);
% %     [bincounts,ind]= histc(concd,binUEdge);
% %     binCent = (binEdge(1:end-1)+binEdge(2:end))./2;
%     
%     %bin and plot concentrations as scatters/bubbles 
%     %x-axis - bins
%     %y-axis - concentrations
% %     concd = conct(:,j);
% %     hsc(j) = scatter(brange,conct(:,j));     
%         
% end



%     set(get(gca(figh),'XLabel'),'String',x_label,...
%                                 'FontName','Arabic Type Setting',...
%                                 'FontSize',12);
%                                     
%     set(get(gca(figh),'YLabel'),'String',y_label,...
%                                 'FontName','Arabic Type Setting',...
%                                 'FontSize',12);   
%     set(get(gca,'XLabel'),'FontWeight','bold');
%     set(get(gca,'YLabel'),'FontWeight','bold');
    

%whitebg(hpss,[0 0 0]);
%function setproperty1(figh,handle,PropStruct)
    
    %figure(handle);
%     if handle == findobj(figh,'type','line')          
%         PropStruct.LineWidth = 1.5;
%         %LpProperty.Color = [.8 .4 0];
%         %PropStruct.LineStyle = 'none';
%         %PropStruct.Marker = 'o';
%         %PropStruct.MarkerSize = 12;
%         %PropStruct.MarkerEdgeColor = [.1 .1 .1];
%         %PropStruct.MarkerFaceColor = [0 .8 0];     
%         
%         
%         set(handle,PropStruct);
%     elseif handle == findall(figh,'type','figure') 
%         PropStruct = struct();
%         %PropStruct.Name = figname;
%         PropStruct.NumberTitle = 'off';
%         %PropStruct.Color = [0.2 0.2 0.2];
%         
%         set(handle,PropStruct);
%         
%     elseif handle == findall(figh,'type','axes')
%         PropStruct.XMinorTick = 'on';
%         %PropStruct.XTick = 0:1000:tmax;
%         PropStruct.YMinorTick = 'on';
%         PropStruct.TickLength = [0.02,0.02];
%         PropStruct.XColor = [.2 .2 .2];
%         PropStruct.YColor = [.2 .2 .2];
%         %PropStruct.XLim = [0 tmax];
%         PropStruct.FontName = 'Courier';
%         PropStruct.FontSize = 12;
%         PropStruct.YGrid = 'on';        
%         %PropStruct.XTickLabel = xticklabel;
%         %PropStruct.XLim = [0 ncols+1];
%         
    
%end
return
function hsubfig = lineplot(j,nplots,hfig,hsubfig,valt,pvalt,tf,y_label1,ColorSpec)
    if hsubfig(j) ~= 0
        hca = findobj(hsubfig(j),'type','axes');
        set(hfig,'CurrentAxes',hca);  
    else%subplot is unassigned
        hsubfig(j) = subplot(2,nplots,j);       
        hca = gca;
        %Make sure more plots can be added at end of loop
        set(hca,'NextPlot','add');     
        set(hsubfig(j),'NextPlot','add');
    end
    
    %Plot only initial and final(sample) concentrations/fluxes
    concd = valt(:,tf);
    pconcd = pvalt(:,tf);   
    X = zeros(length(concd)+length(pconcd),1);        
    initc = concd(1);
    concd(1) = [];
    X(1) = 1;
    X(2:length(pconcd)+1) = 1;
    X(length(pconcd)+2:length(concd)+length(pconcd)+1) = 2;
    hline = line(X,[initc;pconcd;concd]);
    
    %Object Properties
    set(hline,'LineStyle','none',...
              'Marker','o',...
              'MarkerSize',5,...
              'MarkerFaceColor',ColorSpec{j},...
              'MarkerEdgeColor',ColorSpec{j});
    %Axis Properties
    set(get(gca,'YLabel'),'String',y_label1);  
    set(get(gca,'YLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'YLabel'),'FontSize',12); 
    
    xlabel = sprintf('Steady States');
    set(get(gca,'XLabel'),'String',xlabel);  
    set(get(gca,'XLabel'),'FontName','Arabic Type Setting');   
    set(get(gca,'XLabel'),'FontSize',12);
    %Axis limits
    Ydata = get(hline,'YData');
    maxY = max(Ydata);
    minY = min(Ydata);
    [minY,maxY] = fixaxisbounds(minY,maxY);
    
    set(gca,'YLim',[minY maxY]);
    AP.XColor = [.1 .1 .1];
    AP.YColor = [.1 .1 .1];
    AP.XLim = [0 max(X)+1];
    AP.XTick = 0:1:max(X);
    set(hsubfig(j),AP);    
    
return
