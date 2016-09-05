function [flux,vflux] = CKinetics(model,pvec,M,Vind)
[~,nc] = size(M);
allmc = M;
S = model.S;
SI = model.SI;
nrxn = model.nt_rxn;
rev = model.rev;
K = pvec.K;
KIact = pvec.KIact;
KIihb = pvec.KIihb;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,1);
flux = zeros(nrxn,1);

vecmc = repmat(M,1,nrxn);

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



%         sbid = S(:,Vind(irxn))<0;
%         prid = S(:,Vind(irxn))>0;   
%         
%         %remove water
%         sbid(h2o) = 0;
%         prid(h2o) = 0;
%         
%         %remove protons
%         sbid([he hc]) = 0;
%         prid([he hc]) = 0;
%         

%         elseif ~rev(Vind(irxn))                           
%             if all(mc(sbid)>0)
%                 nr_flux =...
%                 kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
%             end
%             if any(sbid) && any(prid)
%                 %Denominator - 1.6
%                 dr_sb = mc(sbid)./K(sbid,Vind(irxn));                             
%                 for j = 1:length(find(sbid))
%                     for si = 2:Sb(j)
%                         dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
%                     end
%                     dr_sb(j) = dr_sb(j) + 1;
%                 end                
%             else
%                 dr_sb = 1;                
%             end
%             dr_pr = 1;
%         end
%       

% non vectorized version for ease of understanding and jacobina
% implementation
for irxn = 1:nrxn
    if ismember(irxn,Vind)
        alls = S(:,irxn);allp = S(:,irxn);
        alls(S(:,irxn)>0) = 0;allp(S(:,irxn)<0) = 0;
        sratio = vecmc(logical(alls),irxn)./K(logical(alls),irxn);
        pratio = vecmc(logical(allp),irxn)./K(logical(allp),irxn);
        thetas = prod(sratio.^-alls(logical(alls)));
        thetap = prod(pratio.^-allp(logical(allp)));
        fwdflx = kfwd(irxn)*thetas;
        revflx = kbkw(irxn)*thetap;
        if ~rev(irxn)
            revflx = 0;
        end
        nrflx = fwdflx-revflx;
        drflx = 1+thetas+thetap;
        
        % partial activation
        if any(SI(:,irxn)>0)
            allac = SI(:,irxn);allac(SI(:,irxn)<0) = 0;
            acratio = vecmc(logical(allac),irxn)./KIact(logical(allac),irxn);
            acflx = prod((0.5+(1-0.5).*acratio./(1+acratio)).^allac(logical(allac)));
        else
            acflx = 1;
        end
        
        % partial inhibition
        if any(SI(:,irxn)<0)
            allib = SI(:,irxn);allib(SI(:,irxn)>0) = 0;
            ibratio = vecmc(logical(allib),irxn)./KIihb(logical(allib),irxn);
            ibflx = prod((0.5+(1-0.5)./(1+ibratio)).^-allib(logical(allib)));
        else
            ibflx = 1;
        end
        nrflx = acflx*ibflx*nrflx;
        
        % specific activation
%         if any(SI(:,irxn)>0)
%             spac = SI(:,irxn);spac(SI(:,irxn)<0) = 0;
%             sparatio = vecmc(logical(spac),irxn)./KIact(logical(spac),irxn);
%             dracflx = sum((1./sparatio).^spac(logical(spac)));
%         else
            dracflx = 0;
%         end
        
        % specific inhibition
%         if any(SI(:,irxn)>0)
%             spib = SI(:,irxn);spib(SI(:,irxn)>0) = 0;
%             spiratio = vecmc(logical(spib),irxn)./KIihb(logical(spib),irxn);
%             dribflx = sum((spiratio).^-spib(logical(spib)));
%         else
            dribflx = 0;
%         end
        drflx = drflx + dracflx + dribflx;
        vflux(irxn) = scale_flux(nrflx/drflx);
        flux(irxn) = Vmax(irxn)*vflux(irxn); 
    end
end

flux = flux(Vind,:);
vflux = vflux(Vind,:);