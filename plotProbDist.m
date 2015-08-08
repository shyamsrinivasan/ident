%Plot probability distribution of fluxes/concentrations
function plotProbDist(model,var,type,varlist)

if nargin < 4
    if strcmpi(type,'flux')
        varlist = model.rxns;        
    elseif strcmpi(type,'mets')
        varlist = model.mets;
    end
end
if strcmpi(type,'flux')
    cmplist = model.rxns;
elseif strcmpi(type,'mets')
    cmplist = model.mets;
end

nvar = length(varlist);
for ivar = 1:nvar    
    if any(strcmpi(varlist{ivar},cmplist))
        tf = strcmpi(varlist{ivar},cmplist);
    %lognormal distribution
        pd = fitdist(var(tf,:)','kernel');
        yvar = pdf(pd,1:.1:100);
        figure
        hline = line(1:.1:100,yvar);
        if strcmpi(type,'mets');
            x_label = sprintf('%s',model.mets{ivar});
        elseif strcmpi(type,'rxns')
            x_label = sprintf('%s',model.rxns{ivar});
        end        
        set(get(gca,'XLabel'),'String',x_label); 
        set(get(gca,'XLabel'),'FontName','CMU Serif');
        set(get(gca,'XLabel'),'FontSize',24);
    end
end

