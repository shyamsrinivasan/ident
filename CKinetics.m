function [flux,vflux] = CKinetics(model,pvec,mc,Vind)
[~,nc] = size(mc);
allmc = mc;
S = model.S;
SI = model.SI;
nrxn = model.nt_rxn;
K = pvec.K;
KIact = pvec.KIact;
KIihb = pvec.KIihb;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,nc);   
flux = zeros(nrxn,nc);

%eliminate consideration for excess cofators
%pi[c],pi[e],h[c],h[e],h2o[c]
% find(strcmpi(model.mets,'pi[e]'))...
%         find(strcmpi(model.mets,'pi[c]'))...
he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
% vmet = [he...
%         find(strcmpi(model.mets,'h[c]'))...        
%         h2o...       
%         find(strcmpi(model.mets,'co2[c]'))];


for ic = 1:nc   
    mc = allmc(:,ic);
    for irxn = 1:length(Vind)
        nr_flux = zeros(1,1);
%         nmet = size(S,1);
        sbid = S(:,Vind(irxn))<0;
        prid = S(:,Vind(irxn))>0;   
        
        %remove water
        sbid(h2o) = 0;
        prid(h2o) = 0;
        
        %remove protons
        sbid([he hc]) = 0;
        prid([he hc]) = 0;
        
        Sb = -S(sbid,Vind(irxn));
        Sp = S(prid,Vind(irxn));

        if model.rev(Vind(irxn))          
            if all(mc(sbid)>0) && all(mc(prid)>0)
              nr_flux = kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb) -...
                        kbkw(Vind(irxn))*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
            elseif all(mc(sbid)>0)
                nr_flux = kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
            elseif all(mc(prid)>0)
                nr_flux = -kbkw(Vind(irxn))*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
            end    
            if any(sbid) && any(prid)
                %Denominator - 1.6
                dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
                for j = 1:length(find(sbid))
                    for si = 2:Sb(j)
                        dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                    end
                end
                %dr_pr
                dr_pr = 1+mc(prid)./K(prid,Vind(irxn));
                for j = 1:length(find(prid))
                    for si = 2:Sp(j)
                        dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
                    end
                end
            else
                dr_sb = 1;
                dr_pr = 1;
            end
        elseif ~model.rev(Vind(irxn))                           
            if all(mc(sbid)>0)
                nr_flux =...
                kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
            end
            if any(sbid) && any(prid)
                %Denominator - 1.6
                dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
                for j = 1:length(find(sbid))
                    for si = 2:Sb(j)
                        dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                    end
                end                  
            else
                dr_sb = 1;                
            end
            dr_pr = 1;
        end
        
        % regulation        
        if any(SI(:,Vind(irxn)))
            % activation
            if any(SI(:,Vind(irxn))>0)
                acid = SI(:,Vind(irxn))>0;
                sac = SI(acid,Vind(irxn));
                nr_flux = nr_flux*prod((mc(acid)./KIact(acid,Vind(irxn))).^sac./...
                          (1+(mc(acid)./KIact(acid,Vind(irxn))).^sac));        
            end
            % inhibition
            if any(SI(:,Vind(irxn))<0)
                ihid = SI(:,Vind(irxn))<0;
                sih = SI(ihid,Vind(irxn));
                nr_flux =...
                nr_flux*prod(1./(1+(mc(ihid)./KIihb(ihid,Vind(irxn))).^sih));
            end
        end  

        if any(sbid) && any(prid)
            dr_flux = prod(dr_sb)+prod(dr_pr)-1;
            vflux(Vind(irxn),ic) = scale_flux(nr_flux/dr_flux);
        else
            vflux(Vind(irxn),ic) = 0;
        end
        flux(Vind(irxn),ic) = Vmax(Vind(irxn))*vflux(Vind(irxn),ic);
    end
end
flux = flux(Vind,:);
vflux = vflux(Vind,:);