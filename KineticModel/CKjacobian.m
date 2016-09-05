function [] = CKjacobian(model,pvec,M,Vind)
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
%             allac = SI(:,irxn);allac(SI(:,irxn)<0) = 0;
%             acratio = vecmc(logical(allac),irxn)./KIact(logical(allac),irxn);
%             dracflx = sum((1./acratio).^allac(logical(allac)));
%         else
            dracflx = 0;
%         end
        
        % specific inhibition
%         if any(SI(:,irxn)>0)
%             allib = SI(:,irxn);allib(SI(:,irxn)>0) = 0;
%             ibratio = vecmc(logical(allib),irxn)./KIihb(logical(allib),irxn);
%             dribflx = sum((ibratio).^-allib(logical(allib)));
%         else
            dribflx = 0;
%         end
        drflx = drflx + dracflx + dribflx;
        vflux(irxn) = scale_flux(nrflx/drflx);
        flux(irxn) = Vmax(irxn)*vflux(irxn); 
        
        % jacobian information
        allTr = zeros(length(find(S(:,irxn))),1);
        allDr = zeros(length(find(S(:,irxn))),1);
        allfr = zeros(length(find(S(:,irxn))),1);
        allDreg = zeros(length(find(S(:,irxn))),1);
        allmet = S(:,irxn)~=0;
        
        zetarxn = fwdflx./revflx;
        
        % numerator dTr/dci,dDr/dci
        metid = find(allmet);
        for im = 1:length(metid)
            if S(metid(im),irxn)<0
                allTr(metid(im)) = (zetarxn*-S(metid(im),irxn)-0)/(zetarxn-1);
                allDr(metid(im)) = thetas*-S(metid(im),irxn)+0*thetap;
            elseif S(metid(im),irxn)>0
                allTr(metid(im)) = (zetarxn*0-S(metid(im),irxn))/(zetarxn-1);
                allDr(metid(im)) = thetas*0+thetap*S(metid(im),irxn);
            end
            allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
            allDr(metid(im)) = allDr(metid(im))/vecmc(metid(im),irxn);
        end
        
        % denominator dDr/dci
        
    end
end

    
% partial derivatives for jacobians
% get all substrates and products
zetar = fwdflx./revflx;
mris = -alls(logical(alls));
mrip = allp(logical(allp));
dnrdc = 1./vecmc(logical(alls)).*(zetar.*mris-mrip)./(zetar-1);

thetas = prod(sratio.^-alls(logical(alls)),2);
thetap = prod(pratio.^allp(logical(allp)),2);
ddrdc = 1/ci*(mris*thetas - mrip*thetap)/(thetas+thetap+1);
% dv/dc = vr*(1/fr*dfr/dc + 1/num*dnum/dc+1/dem*ddem/dc)