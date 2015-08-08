function ensb = resample_ensemble(ensb,model,MC)
nmodels = length(fieldnames(ensb));

%glcpts
Mglc = strcmpi('glc[e]',model.mets);
Vglc = find(model.S(Mglc,:)<0);

try
    Vind = [model.Vind;Vglc];
catch
    Vind = [model.Vind Vglc];
end
nrxn = length(Vind);

%Order of magnitude to determine whether to resample Km or not
ofm = 1;
rfmlb = -1;
rfmub = 1;
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
        %Regulatory parameters
        if any(model.SI(:,Vind(irxn))~=0)                
            rKAind = model.SI(:,Vind(irxn))>0;%activation
            rpKAjind = find(set.KAind(rKAind,Vind(irxn))~=1);
            rpKA = set.KIact(rKAind,Vind(irxn));
            rpAC = MC(rKAind);                   
            %order of magnitude for KIact based on C
            %concentration +/- 100
            if any(round(log10(rpKA)-log10(rpAC))>rfmub) ||...
               any(round(log10(rpKA)-log10(rpAC))<rfmlb)
                indlb = round(log10(rpKA)-log10(rpAC))<rfmlb;
                indub = round(log10(rpKA)-log10(rpAC))>rfmub;
                ind = setdiff([find(indlb) find(indub)],rpKAjind);
                while any(round(log10(rpKA(ind))-log10(rpAC(ind)))>rfmub) ||...
                      any(round(log10(rpKA(ind))-log10(rpAC(ind)))<rfmlb)                
                    %new sample
                    reg_sat = betarnd(1.5,4.5,length(find(ind)),1);
                    %new ratio
                    act_ratio = reg_sat./(1-reg_sat);
                    rpKA(ind) = rpAC(ind)./act_ratio;
                end                    
            end
            set.KIact(rKAind,Vind(irxn)) = rpKA;
            rKIind = model.SI(:,Vind(irxn))<0;%inhibition
            rpKIjind = find(set.KIind(rKIind,Vind(irxn))~=1);
            rpKI = set.KIihb(rKIind,Vind(irxn));
            rpIC = MC(rKIind); 
            %order of magnitude for KIihb based on C
            %concentration +/- 100            
            if any(round(log10(rpKI)-log10(rpIC))>rfmub) ||...
               any(round(log10(rpKI)-log10(rpIC))<rfmlb)
                indlb = round(log10(rpKI)-log10(rpIC))<rfmlb;
                indub = round(log10(rpKI)-log10(rpIC))>rfmub;
                ind = setdiff([find(indlb) find(indub)],rpKIjind); 
                while any(round(log10(rpKI(ind))-log10(rpIC(ind)))>rfmub) ||...
                      any(round(log10(rpKI(ind))-log10(rpIC(ind)))<rfmlb)      
                    %new sample
                    reg_sat = betarnd(1.5,4.5,length(find(ind)),1);
                    %new ratio
                    neg_ratio = reg_sat./(1-reg_sat);
                    rpKI(ind) = rpIC(ind)./neg_ratio;
                end
            end
            set.KIihb(rKIind,Vind(irxn)) = rpKI;
        end
    end
    %calculate new Vmax
    [~,vflux] = ConvinienceKinetics(model,set,MC,model.bmrxn,Vind);
    oldVmax = set.Vmax;
    set.Vmax(Vind) = model.Vss(Vind)./vflux(Vind);    
    set.Vmax(setdiff(1:model.nt_rxn,Vind)) = 1; 
    if any(~isnan(oldVmax))
        set.Vmax(~isnan(oldVmax)) = oldVmax(~isnan(oldVmax));
    end
    ensb.(mname) = set;
    ensb.(mname) = rmfield(ensb.(mname),{'Kind','KAind','KIind'});
end


%Select models from ensemble
for is = 1:nmodels
     mname = sprintf('model%d',is);  
    if any(ensb.(mname).Vmax < 0)
        fprintf('%s not suitable for simulation\n',mname);
    end
end

%Select models from ensemble
for is = 1:nmodels
     mname = sprintf('model%d',is);  
    if any(ensb.(mname).Vmax < 0)
        fprintf('%s not suitable for simulation\n',mname);
    end
end

return