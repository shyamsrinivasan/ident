function flux = Kotte2014Ckinetics_flux(allmc,model,pvec)
tfenz = strcmpi(model.mets,'enz[c]');
flux = CKinetics(model,pvec,allmc,[1 2 3 4]);
flux(strcmpi(model.rxns,'ACpts')) = allmc(tfenz)*flux(strcmpi(model.rxns,'ACpts'));
flux = Kotte_glycolysisflux(allmc,pvec,flux,model);    