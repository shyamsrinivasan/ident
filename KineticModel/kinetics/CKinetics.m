function [flux,vflux,DVX] = CKinetics(model,pvec,M,Vind,getjac)
if nargin<5
    getjac=0;
end
[~,nc] = size(M);
S = model.S;
SI = model.SI;
nrxn = model.nt_rxn;
rev = model.rev;
remid = model.remid;
Mext = model.Mext;
rho = model.rho;
biomass = strcmpi(model.mets,'biomass[e]');

% convert units mmole/Lc -> mmole/Lcw
M(Mext) = M(Mext).*rho/M(biomass);

K = pvec.K;
KIact = pvec.KIact;
KIihb = pvec.KIihb;
kfwd = pvec.kfwd;
kbkw = pvec.krev;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,nc);
flux = zeros(nrxn,nc);
DVX = zeros(size(M,1),nrxn);

% vflux = cons(vflux,M);
% flux = cons(flux,M);

% vecmc = repmat(M,1,nrxn);

% eliminate consideration for excess cofators
% pi[c],pi[e],h[c],h[e],h2o[c]
% he = find(strcmpi(model.mets,'h[e]'));
% hc = find(strcmpi(model.mets,'h[c]'));
% h2o = find(strcmpi(model.mets,'h2o[c]'));
% pi = find(strcmpi(model.mets,'pi[c]'));

% find(strcmpi(model.mets,'pi[e]'))...
%         find(strcmpi(model.mets,'pi[c]'))...

%                 %Denominator - 1.6
%                 dr_sb = mc(sbid)./K(sbid,Vind(irxn));                             
%                 for j = 1:length(find(sbid))
%                     for si = 2:Sb(j)
%                         dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
%                     end
%                     dr_sb(j) = dr_sb(j) + 1;
%                 end                


% non vectorized version for ease of understanding and jacobina
% implementation
for irxn = 1:nrxn
    if ismember(irxn,Vind)
        alls = S(:,irxn);allp = S(:,irxn);
        alls(S(:,irxn)>0) = 0;allp(S(:,irxn)<0) = 0;
        alls(remid,:) = 0;allp(remid,:) = 0;
        if nc>1
            sratio = M(logical(alls),:)./...
                     repmat(K(logical(alls),irxn),1,nc);
            pratio = M(logical(allp),:)./...
                     repmat(K(logical(allp),irxn),1,nc);
            thetas = prod(sratio.^...
                     repmat(-alls(logical(alls)),1,nc),1);
            thetap = prod(pratio.^...
                     repmat(allp(logical(allp)),1,nc),1);
        else
            sratio = M(logical(alls),:)./K(logical(alls),irxn);
            pratio = M(logical(allp),:)./K(logical(allp),irxn);
            thetas = prod(sratio.^-alls(logical(alls)),1);
            thetap = prod(pratio.^allp(logical(allp)),1);
        end    
        fwdflx = kfwd(irxn).*thetas;
        revflx = kbkw(irxn).*thetap;
        if ~rev(irxn)
%             revflx = zeros(1,length(revflx));
        end
        nrflx = fwdflx-revflx;
        drflx = 1+thetas+thetap;
        
        % partial activation
        rhoA = 0.5;
        if any(SI(:,irxn)>0)
            allac = SI(:,irxn);allac(SI(:,irxn)<0) = 0;
            if nc>1
                acratio = M(logical(allac),:)./...
                          repmat(KIact(logical(allac),irxn),1,nc);
                acflx = prod((rhoA+(1-rhoA).*acratio./(1+acratio)).^...
                        repmat(allac(logical(allac)),1,nc),1);
            else
                acratio = M(logical(allac),:)./KIact(logical(allac),irxn);
                acflx = prod((rhoA+(1-rhoA).*acratio./(1+acratio)).^...
                        allac(logical(allac)),1);
            end            
        else
            allac = [];
            acratio = [];
            acflx = ones(1,length(nrflx));
        end
        
        % partial inhibition
        rhoI = 0.5;
        if any(SI(:,irxn)<0)
            allib = SI(:,irxn);allib(SI(:,irxn)>0) = 0;
            if nc>1
                ibratio = M(logical(allib),:)./...
                          repmat(KIihb(logical(allib),irxn),1,nc);
                ibflx = prod((rhoI+(1-rhoI)./(1+ibratio)).^...
                        repmat(-allib(logical(allib)),1,nc),1);
            else
                ibratio = M(logical(allib),:)./KIihb(logical(allib),irxn);
                ibflx = prod((rhoI+(1-rhoI)./(1+ibratio)).^...
                        -allib(logical(allib)),1);
            end
        else
            allib = [];
            ibratio = [];
            ibflx = ones(1,length(nrflx));
        end
        nrflx = acflx.*ibflx.*nrflx;
        
        % specific activation
%         if any(SI(:,irxn)>0)
%             spac = SI(:,irxn);spac(SI(:,irxn)<0) = 0;
%             sparatio = vecmc(logical(spac),irxn)./KIact(logical(spac),irxn);
%             dracflx = sum((1./sparatio).^spac(logical(spac)));
%         else
            dracflx = zeros(1,length(drflx));
%         end
        
        % specific inhibition
%         if any(SI(:,irxn)>0)
%             spib = SI(:,irxn);spib(SI(:,irxn)>0) = 0;
%             spiratio = vecmc(logical(spib),irxn)./KIihb(logical(spib),irxn);
%             dribflx = sum((spiratio).^-spib(logical(spib)));
%         else
            dribflx = zeros(1,length(drflx));
%         end
        drflx = drflx + dracflx + dribflx;
        vflux(irxn,:) = scale_flux(nrflx./drflx);
        flux(irxn,:) = Vmax(irxn).*vflux(irxn,:); 
        
        % get jacobian information
        if getjac
            DVX(:,irxn) = getCKjacobian(irxn,flux,S,SI,rev,sratio,pratio,thetas,thetap,...
                                        nrflx,drflx,fwdflx,revflx,...
                                        allac,acratio,rhoA,allib,ibratio,rhoI);
        end
    end
end

if getjac
    DVX = DVX(:,Vind);
end

flux = flux(Vind,:);
vflux = vflux(Vind,:);