function flux  = changefluxunits(flux,M,model,convert)
[nf,nc] = size(flux);
biomass = strcmpi(model.mets,'biomass[e]');
rho = model.rho;
if nc>1
    switch convert
        case 1 % mmole/Lcw.h -> mmole/Lc.h     
            flux = flux.*repmat(M(biomass,:),nf,1)./repmat(rho,nf,nc);
        case 2 % mmole/gDCW.h -> mmole/Lcw.h
            flux = flux.*repmat(rho,nf,nc);
        case 3 % mmole/gDCW.h -> mmole/Lc.h
            flux = flux.*repmat(M(biomass,:),nf,1);
        case -1 % mmole/Lc.h -> mmole/Lcw.h
            flux = flux.*repmat(rho,nf,nc)./repmat(M(biomass,:),nf,1);
        case -2 % mmole/Lcw.h -> mmole/gDCW.h
            flux = flux./repmat(rho,nf,nc);
        case -3 % mmole/Lc.h -> mmole/gDCW.h
            flux = flux./repmat(M(biomass,:),nf,1);
    end
else
    switch convert
        case 1 % mmole/Lcw.h -> mmole/Lc.h     
            flux = flux.*M(biomass)./rho;
        case 2 % mmole/gDCW.h -> mmole/Lcw.h
            flux = flux.*rho;
        case 3 % mmole/gDCW.h -> mmole/Lc.h
            flux = flux.*M(biomass);
        case -1 % mmole/Lc.h -> mmole/Lcw.h
            flux = flux.*rho./M(biomass);
        case -2 % mmole/Lcw.h -> mmole/gDCW.h
            flux = flux./rho;
        case -3 % mmole/Lc.h -> mmole/gDCW.h
            flux = flux./M(biomass);
    end
end