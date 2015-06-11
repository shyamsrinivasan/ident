function dYdt = changeC(t,EC,Mrel,dYdt,model,ng,gr_flux)
if t > 360000 && t < 720000
    mt_ind = find(strcmpi('P',model.Regulators))+ng(1);   
    mx_ind = strcmpi('P',model.Metabolites);
    dYdt(mt_ind) = .0005/(1+exp(-1*(t-360000)))- 0.1*Mrel(mx_ind)*gr_flux;
%     dYdt(mt_ind) = - 0.1*Mrel(mx_ind)*gr_flux;
end
return