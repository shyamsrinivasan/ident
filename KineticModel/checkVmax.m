function flag = checkVmax()
pvector.K = Km;
pvector.KIact = KIact;
pvector.KIihb = KIihb;
pvector.kact_ratio = 
pvector.kcat_fwd
pvector.Vmax = 1;
[~,vflux] = ConvinienceKinetics(model,pvector,MC,model.bmrxn,irxn);