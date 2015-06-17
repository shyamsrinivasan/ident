function growth = biomass_flux(model,Y,dXdt,flux)
%Metabolites consumed for biomass
Mb = model.S(:,model.bmrxn)<0;
ATP = strcmpi('atp[c]',model.Metabolites);
Mb = setdiff(find(Mb),find(ATP));

Sj = -model.S(Mb,model.bmrxn);
Mj = model.MolWt(Mb);%Need to get this done

dXdt(Mb) = model.S(Mb,:)*flux;
if all(Y(Mb))
    growth = sum((Sj.*Mj.*1e-3.*Y(Mb)))/(1+sum(Sj.*Mj.*1e-3.*Y(Mb)));
else
    growth = 0;
end
    


