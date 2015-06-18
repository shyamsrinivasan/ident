function ensb = resample_ensemble(ensb,model,MC)
nmodels = length(fieldnames(ensb));

%glcpts
Mglc = strcmpi('glc[e]',model.Metabolites);
Vglc = find(model.S(Mglc,:)<0);

try
    Vind = [model.Vind;Vglc];
catch
    Vind = [model.Vind Vglc];
end
nrxn = length(Vind);

%Order of magnitude to determine whether to resample Km or not
ofm = 1;
for isample = 1:nmodels
    mname = sprintf('model%d',isample);
    set = ensb.(mname);    
    for irxn = 1:nrxn        
        rpind = model.S(:,Vind(irxn))~=0;        
        rpKjnd = find(set.Kind(rpind,Vind(irxn))~=1);
        rpK = set.K(rpind,Vind(irxn));
        rpC = MC(rpind);
        %order of magnitude for Km based on C
        %concentration +/- 100
        if any(round(log10(rpK)-log10(rpC))>ofm)
            ind = round(log10(rpK)-log10(rpC))>ofm;
            ind = setdiff(find(ind),rpKjnd);
            while any(round(log10(rpK(ind))-log10(rpC(ind)))>1)
                %resample saturation rate
                met_sat = betarnd(1.5,4.5,length(find(ind)),1);
                sub_ratio = met_sat./(1-met_sat); 
                rpK(ind) = rpC(ind)./sub_ratio;
            end
        end
        set.K(rpind,Vind(irxn)) = rpK;        
    end
    %calculate new Vmax
    [~,vflux] = ConvinienceKinetics(model,set,MC,model.bmrxn,Vind);
    set.Vmax(Vind) = model.Vss(Vind)./vflux(Vind);    
    set.Vmax(setdiff(1:model.nt_rxn,Vind)) = 1; 
    ensb.(mname) = set;
end

%Select models from ensemble
for is = 1:nmodels
     mname = sprintf('model%d',is);  
    if any(ensb.(mname).Vmax < 0)
        fprintf('%s not suitable for simulation\n',mname);
    end
end

return