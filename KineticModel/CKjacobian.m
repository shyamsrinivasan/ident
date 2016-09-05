function [DVX] = CKjacobian(model,pvec,M,Vind)
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

vflux = zeros(nrxn,nc);
flux = zeros(nrxn,nc);
DVX = zeros(length(M),nrxn);

vecmc = repmat(M,1,nrxn);

% vectorization does not work as well wrt reactions
% switching back to using for loops for rxn
for irxn = 1:nrxn
    if ismember(irxn,Vind)
        alls = S(:,irxn);allp = S(:,irxn);
        alls(S(:,irxn)>0) = 0;allp(S(:,irxn)<0) = 0;
        sratio = vecmc(logical(alls),irxn)./K(logical(alls),irxn);
        pratio = vecmc(logical(allp),irxn)./K(logical(allp),irxn);
        thetas = prod(sratio.^-alls(logical(alls)));
        thetap = prod(pratio.^allp(logical(allp)));
        fwdflx = kfwd(irxn)*thetas;
        revflx = kbkw(irxn)*thetap;
        if ~rev(irxn)
            revflx = 0;
        end
        nrflx = fwdflx-revflx;
        drflx = 1+thetas+thetap;
        
        % partial activation
        rhoA = 0.5;
        if any(SI(:,irxn)>0)
            allac = SI(:,irxn);allac(SI(:,irxn)<0) = 0;
            acratio = vecmc(logical(allac),irxn)./KIact(logical(allac),irxn);
            acflx = prod((rhoA+(1-rhoA).*acratio./(1+acratio)).^allac(logical(allac)));
        else
            acflx = 1;
        end
        
        % partial inhibition
        rhoI = 0.5;
        if any(SI(:,irxn)<0)
            allib = SI(:,irxn);allib(SI(:,irxn)>0) = 0;
            ibratio = vecmc(logical(allib),irxn)./KIihb(logical(allib),irxn);
            ibflx = prod((rhoI+(1-rhoI)./(1+ibratio)).^-allib(logical(allib)));
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
        
        % jacobian information
        allTr = zeros(length(S(:,irxn)),1);
        allDr = zeros(length(S(:,irxn)),1);
        allfr = zeros(length(SI(:,irxn)),1);
        allDreg = zeros(length(SI(:,irxn)),1);
        allvr = zeros(length(S(:,irxn)),1);
        allmet = S(:,irxn)~=0;
        metid = find(allmet);
        
        % numerator dTr/dci
        if rev(irxn) && all(pratio>0)
            zetarxn = fwdflx./revflx;
            for im = 1:length(metid)
                if S(metid(im),irxn)<0
                    allTr(metid(im)) = (zetarxn*-S(metid(im),irxn)-0)/(zetarxn-1);                    
                elseif S(metid(im),irxn)>0
                    allTr(metid(im)) = (zetarxn*0-S(metid(im),irxn))/(zetarxn-1);                    
                end
                allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);                
            end
        elseif ~rev(irxn)
            for im = 1:length(metid)
                if S(metid(im),irxn)<0
                    allTr(metid(im)) = -S(metid(im),irxn);
                end
                allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
            end
        end
        
        % denominator dDr/dci        
        for im = 1:length(metid)
            if S(metid(im),irxn)<0                
                allDr(metid(im)) = thetas*-S(metid(im),irxn)+0*thetap;
            elseif S(metid(im),irxn)>0                
                allDr(metid(im)) = thetas*0+thetap*S(metid(im),irxn);
            end            
            allDr(metid(im)) = allDr(metid(im))/vecmc(metid(im),irxn);
        end
        
        % denominator specific regulation dDreg/dci
        % specific activation   
        if any(SI(:,irxn)>0)
%             allDreg(logical(spac)) = spac(logical(spac)).*...
%                                      (1./spiratio).^spib(logical(spac))./...
%                                      vecmc(logical(spac),irxn);
        end
        
        % specific inhibition
        if any(SI(:,irxn)<0)
%             allDreg(logical(spib)) = spib(logical(spib)).*...
%                                      (spiratio).^spib(logical(spib))./...
%                                      vecmc(logical(spib),irxn);
        end
        
        allDr = allDr./drflx;
        allDreg = allDreg./drflx;
        
        % partial/essential activation dacflx/dci        
        if any(SI(:,irxn)>0)
            alphaA = 1./(1+acratio);
            betaA = acratio.*alphaA;
            allfr(logical(allac)) =...
            (1-rhoA)./(1-rhoA+rhoA.*betaA.^(-allac(logical(allac)))).*...
            allac(logical(allac)).*alphaA;            
        end
        
        % partial/essential inhibition dibflx/dci        
        if any(SI(:,irxn)<0)
            alphaI = 1./(1+ibratio);
            betaI = ibratio.*alphaI;
            allfr(logical(allib)) = allfr(logical(allib))...
            -(1-rhoI)./(1-rhoI+rhoI.*alphaI.^(-allib(logical(allib)))).*...
            allib(logical(allib)).*betaI;
        end
        
        % collect common ids between allac and allid
        if any(SI(:,irxn)>0) && any(SI(:,irxn)<0)
            cmnreg = union(find(allac),find(allib));
            % remove common ids from allac and allib
            allac = allac(~ismember(allac,cmnreg));
            allib = allib(~ismember(allib,cmnreg));
            allfr(cmnreg) = allfr(cmnreg)./vecmc(cmnreg,irxn);
        end
        if any(SI(:,irxn)>0)
            if ~isempty(allac)
                allfr(logical(allac)) = allfr(logical(allac))./vecmc(logical(allac),irxn);
            end
        end
        if any(SI(:,irxn)<0)
            if ~isempty(allib)
                allfr(logical(allib)) = allfr(logical(allib))./vecmc(logical(allib),irxn);
            end  
        end
        
        % complete differential w.r.t flux(irxn) dvr/dci
        allvr = flux(irxn).*(allfr + allTr - allDr - allDreg);
        DVX(:,irxn) = allvr;        
    end
end
DVX = DVX(:,Vind);

    
% partial derivatives for jacobians
% get all substrates and products

% dv/dc = vr*(1/fr*dfr/dc + 1/num*dnum/dc+1/dem*ddem/dc)