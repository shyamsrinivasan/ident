%PLot kinetic parameters 
function plotParameters(model,ensb)
nmodels = length(fieldnames(ensb));
K = zeros(length(model.mets),model.nt_rxn,nmodels);

for imodel = 1:nmodels
    mname = sprintf('model%d',imodel);     
    for irxn = 1:model.nt_rxn
        K(:,irxn,imodel) = full(ensb.(mname).K(:,irxn));        
    end
end
Vind = model.Vind;
for irxn = 1:length(Vind)
    [subs,~] = find(model.S(:,Vind(irxn)));    
    for imet = 1:length(subs)
        Kval = reshape(K(subs(imet),Vind(irxn),:),...
                       size(K(subs(imet),Vind(irxn),:),2),[]);
        minK = min(Kval);
        maxK = max(Kval);
        pd = fitdist(Kval','kernel');
        yvar = pdf(pd,minK:.1:maxK);
        figure
        hline = line(minK:.1:maxK,yvar);   
        x_label = sprintf('%s, Flux:%s',model.mets{subs(imet)},model.rxns{Vind(irxn)});
        set(get(gca,'XLabel'),'String',x_label); 
        set(get(gca,'XLabel'),'FontName','CMU Serif');
        set(get(gca,'XLabel'),'FontSize',24);
    end
end


