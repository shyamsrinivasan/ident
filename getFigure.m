function getFigure(conc,flux,model)
set(0,'DefaultTextFontname', 'CMU Serif')
set(0,'DefaultAxesFontName', 'CMU Serif')
set(0,'defaulttextinterpreter','latex')
nmet = size(conc,2);
for im = 1:nmet    
    for ij = 1:nmet
        if im ~= ij && ij > im
            figure
            plot(conc(:,im),conc(:,ij),'LineStyle','none',...
                                       'Marker','o',...
                                       'MarkerFaceColor','none',...
                                       'MarkerEdgeColor',[0 0 0]);
            ylabel = sprintf('%s mM',model.mets{ij});
            xlabel = sprintf('%s mM',model.mets{im});
            set(get(gca,'YLabel'),'String',ylabel);  
            set(get(gca,'XLabel'),'String',xlabel); 
        end
    end                           
end
close all
nrxn = size(flux,2);
for im = 1:nrxn  
    for ij = 1:nrxn
        if im ~= ij && ij > im
            if any(strcmpi('Pout',model.rxns{ij})) ||...
               any(strcmpi('Pout',model.rxns{im})) ||...
               any(strcmpi('v6',model.rxns{ij})) ||...
               any(strcmpi('v6',model.rxns{im})) ||...
               any(strcmpi('v8',model.rxns{ij})) ||...
               any(strcmpi('v8',model.rxns{im})) 
                figure
                plot(flux(:,im),flux(:,ij),'LineStyle','none',...
                                           'Marker','o',...
                                           'MarkerFaceColor','none',...
                                           'MarkerEdgeColor',[0 0 0]);
                ylabel = sprintf('%s mmole/gDCW.h',model.rxns{ij});
                xlabel = sprintf('%s mmole/gDCW.h',model.rxns{im});
                set(get(gca,'YLabel'),'String',ylabel);  
                set(get(gca,'XLabel'),'String',xlabel); 
            end
        end
    end                           
end
close all