function DVX = TKjacobian(model,pvec,M,Vex)
S = model.S;
nrxn = model.nt_rxn;
rev = model.rev;
K = pvec.K;
kfwd = pvec.kfwd;
kbkw = pvec.krev;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,1);
flux = zeros(nrxn,1);
DVX = zeros(length(M),nrxn);

vecmc = repmat(M,1,nrxn);

pie = strcmpi(model.mets,'pi[e]');

for irxn = 1:nrxn
    if ismember(irxn,Vex)
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
        % numerator of kinetics
        nrflx = fwdflx-revflx;
        drflx = 1+thetas+thetap;
        vflux(irxn) = scale_flux(nrflx/drflx);
        
        % fluxes for O2t and PIt2r from Vex - overwrite nrflx and drflx
        if irxn == find(strcmpi(model.rxns,'O2t'))
        end
        if irxn == find(strcmpi(model.rxns,'PIt2r'))
            Kapie = 0.89; % mM
            vflux(irxn) =...
            scale_flux(fwdflx/drflx*1/(1+Kapie/vecmc(pie,irxn)));
        end
        flux(irxn) = Vmax(irxn)*vflux(irxn); 
        
        % get jacobian information
        allTr = zeros(length(S(:,irxn)),1);
        allDr = zeros(length(S(:,irxn)),1);
        allvr = zeros(length(S(:,irxn)),1);
        allmet = S(:,irxn)~=0;
        metid = find(allmet);
        
        % numerator dTr/dci
        if rev(irxn) && all(pratio>0)
            zetarxn = fwdflx./revflx;
            for im = 1:length(metid)
                if S(metid(im),irxn)<0
                    if ~isinf(zetarxn)
                        allTr(metid(im)) = (zetarxn*-S(metid(im),irxn)-0)/(zetarxn-1); 
                    else
                        allTr(metid(im)) = 1; 
                    end
                elseif S(metid(im),irxn)>0
                    if ~isinf(zetarxn)
                        allTr(metid(im)) = (zetarxn*0-S(metid(im),irxn))/(zetarxn-1);
                    else
                        allTr(metid(im)) = 1;
                    end
                end
%                 if vecmc(metid(im),irxn)>0
%                     allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
%                 else
%                     allTr(metid(im)) = 0;
%                 end
            end
        elseif ~rev(irxn)
            for im = 1:length(metid)
                if S(metid(im),irxn)<0
                    allTr(metid(im)) = -S(metid(im),irxn);
                end
%                 if vecmc(metid(im),irxn)>0
%                     allTr(metid(im)) = allTr(metid(im))/vecmc(metid(im),irxn);
%                 else
%                     allTr(metid(im)) = 0;
%                 end
            end
        end
        
        % denominator dDr/dci        
        for im = 1:length(metid)
            if S(metid(im),irxn)<0                
                allDr(metid(im)) = thetas*-S(metid(im),irxn)+0*thetap;
            elseif S(metid(im),irxn)>0                
                allDr(metid(im)) = thetas*0+thetap*S(metid(im),irxn);
            end         
            if irxn == find(strcmpi(model.rxns,'PIt2r'))
                if metid(im) == find(pie)
                    piratio = Kapie/vecmc(metid(im),irxn);
                    allDr(metid(im)) = allDr(metid(im))/drflx-...
                    (piratio/vecmc(metid(im),irxn))/(1+piratio);
                end
            else
                allDr(metid(im)) = allDr(metid(im))/drflx;
            end
%             if vecmc(metid(im),irxn)>0
%                 allDr(metid(im)) = allDr(metid(im))/vecmc(metid(im),irxn);
%             else
%                 allDr(metid(im)) =  0;
%             end
        end             
        
        % complete differential w.r.t flux(irxn) dvr/dci
        allvr = flux(irxn).*(allTr - allDr);
        DVX(:,irxn) = allvr;    
    end
end
DVX = DVX(:,Vex);
