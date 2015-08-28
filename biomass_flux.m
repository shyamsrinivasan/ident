function mu = biomass_flux(model,met,dXdt,flux)
%Metabolites consumed for biomass
% Mb = model.S(:,model.bmrxn)<0;
% ATP = strcmpi('atp[c]',model.mets);
% Mb = setdiff(find(Mb),find(ATP));

% dXdt = zeros(length(Mb),1);
% Sj = -model.S(Mb,model.bmrxn);
% Mj = model.MolWt(Mb);%Need to get this done
% YMb = Y(Mb);
% 
% dXdt = model.S(Mb,:)*flux-flux(model.bmrxn)*Y(Mb);
% if all(Y(Mb)>0)
% %     growth = sum((Sj.*Mj.*1e-3.*dXdt(Mb)))/(1+sum(Sj.*Mj.*1e-3.*Y(Mb)));
% %     growth = sum((Sj.*Mj.*1e-3.*Y(Mb)))/(1+sum(Sj.*Mj.*1e-3.*Y(Mb)));
% %set thresholds on Metabolites
%     minflux = min(dXdt);
%     minID = find(dXdt(dXdt==minflux));
%     if dXdt(minID) >=0
%         growth = dXdt(minID)/Sj(minID)/(1+YMb(minID)/Sj(minID));
%     else
%         growth = 0;
%     end
% else
%     growth = 0;
% end
bm = model.S(:,model.bmrxn)<0;
Sj = -model.S(bm,model.bmrxn);
Vnobm = setdiff(1:model.nt_rxn,model.bmrxn);
netflux = zeros(length(Vnobm),1);
netflux(bm) = model.S(bm,Vnobm)*flux(Vnobm);
if all(met(bm)>0)
    vsynth = netflux(bm)./(Sj.*(1+met(bm)./Sj));
    vmin = min(vsynth);
    if vmin >= 0
        mu = vmin;
    else
        mu = -vmin;
    end    
else
    mu = 0;
end
                                                                                                                                                                                        
    


