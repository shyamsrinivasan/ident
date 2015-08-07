function growth = biomass(model,Ynew,Yold,tnew,told,olflux,newflux)

Mbio = strcmpi('biomass',model.mets);

%Metabolites consumed for biomass
Mb = model.S(:,model.bmrxn)<0;
ATP = strcmpi('atp[c]',model.mets);
Mb = setdiff(find(Mb),find(ATP));

Sj = -model.S(Mb,model.bmrxn);
Mj = model.MolWt(Mb);%Need to get this done

if ~isempty(Yold) && ~isempty(told) && ~isempty(told)
    if all(Yold(Mb))
        X0p = sum((Sj.*Mj.*1e-3.*Yold(Mb)));%/(1+sum(Sj.*Mj.*1e-3.*Yold(Mb)));
    else
        X0p = 0;
    end

    if all(Ynew(Mb))
        X1p = sum((Sj.*Mj.*1e-3.*Ynew(Mb)));%/(1+sum(Sj.*Mj.*1e-3.*Ynew(Mb)));
    else
        X1p = 0;
    end
    if tnew > told 
        if X1p > X0p
            growth = (X1p - X0p)/((tnew-told)*X1p);
        else
            growth = 0;
        end
    else
        growth = olflux(model.bmrxn);
    end
else    
    growth = 0;
end



