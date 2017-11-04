function flux = Kotte2014Ckinetics_flux(allmc,model,pvec)
tfenz = strcmpi(model.mets,'enz[c]');
flux([1;3;4]) = CKinetics(model,pvec,allmc,[1 3 4]);
flux(strcmpi(model.rxns,'ACpts')) = allmc(tfenz)*flux(strcmpi(model.rxns,'ACpts'));
flux = Kotte_glycolysisflux(allmc,pvec,flux,model);    