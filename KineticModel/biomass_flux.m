function growth = biomass_flux(model,dXdt,flux)
%Metabolites consumed for biomass
Mb = model.S(:,model.bmrxn)<0;
ATP = strcmpi('atp[c]',model.Metabolites);
Mb = setdiff(find(Mb),find(ATP));

Sj = model.S(Mb,model.bmrxn);
Mj = model.MolWt(Mb,model.bmrxn);%Need to get this done

dXdt(Mb) = model.S(Mb,:)*flux;
growth = sum((Sj.*Mj.*dXdt(Mb)))/(1+sum(Sj.*Mj));


