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

vflux = zeros(nrxn,nc);   
vflux = cons(vflux,M);
flux = zeros(nrxn,nc);
flux = cons(flux,M);

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


% for ic = 1:nc   
%     mc = allmc(:,ic);
%     for irxn = 1:length(Vind)
%         nr_flux = zeros(1,1);
% %         nmet = size(S,1);
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
%         Sb = -S(sbid,Vind(irxn));
%         Sp = S(prid,Vind(irxn));
% 
%         if rev(Vind(irxn))          
%             if all(mc(sbid)>0) && all(mc(prid)>0)
%               nr_flux = kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb) -...
%                         kbkw(Vind(irxn))*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
%             elseif all(mc(sbid)>0)
%                 nr_flux = kfwd(Vind(irxn))*prod((mc(sbid)./K(sbid,Vind(irxn))).^Sb);
%             elseif all(mc(prid)>0)
%                 nr_flux = -kbkw(Vind(irxn))*prod((mc(prid)./K(prid,Vind(irxn))).^Sp);
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
%                 %dr_pr
%                 dr_pr = mc(prid)./K(prid,Vind(irxn));                
%                 for j = 1:length(find(prid))
%                     for si = 2:Sp(j)
%                         dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
%                     end
%                     dr_pr(j) = dr_pr(j) + 1;
%                 end
%             else
%                 dr_sb = 1;
%                 dr_pr = 1;
%             end
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
%         % regulation        
%         if any(SI(:,Vind(irxn)))
%             % activation
%             if any(SI(:,Vind(irxn))>0)
%                 acid = SI(:,Vind(irxn))>0;
%                 sac = SI(acid,Vind(irxn));
%                 acflx = prod((0.5+(1-0.5).*mc(acid)./KIact(acid,Vind(irxn))/...
%                         (1+mc(acid)./KIact(acid,Vind(irxn)))).^sac);
%                 nr_flux = acflx*nr_flux;
% %                 nr_flux = nr_flux*...
% %                           prod(1 + (mc(acid)./KIact(acid,Vind(irxn))).^sac);                      
%             end
%             % inhibition
%             if any(SI(:,Vind(irxn))<0)
%                 ihid = SI(:,Vind(irxn))<0;
%                 sih = SI(ihid,Vind(irxn));
%                 ibflx = prod((0.5+(1-0.5).*1/...
%                         (1+mc(acid)./KIact(acid,Vind(irxn)))).^sac);
%                 nr_flux = ibflx*nr_flux;
% %                 nr_flux =...
% %                 nr_flux*prod(1./(1+(mc(ihid)./KIihb(ihid,Vind(irxn))).^sih));                
%             end
%         end  
% 
%         if any(sbid) && any(prid)
%             dr_flux = prod(dr_sb)+prod(dr_pr)-1;
%             vflux(Vind(irxn),ic) = scale_flux(nr_flux/dr_flux);
%         else
%             vflux(Vind(irxn),ic) = 0;
%         end
%         flux(Vind(irxn),ic) = Vmax(Vind(irxn))*vflux(Vind(irxn),ic);
%     end
% end

% vectorized version of CKinetics for all rxns in Vind - under testing
for ic = 1:nc
    M = allmc(:,ic);
    vecmc = repmat(M,1,length(Vind));
    acflx = ones(model.nt_rxn,1);
    ibflx = ones(model.nt_rxn,1);
    sprod = ones(length(Vind),1);
    pprod = ones(length(Vind),1);
%     acdr = zeros(model.nt_rxn,1);
%     ibdr = zeros(model.nt_rxn,1);
    alls = S(:,Vind);allp = S(:,Vind);
    alls(S(:,Vind)>0) = 0;allp(S(:,Vind)<0) = 0;
    allK = K(:,Vind);
    sratio = vecmc(logical(alls))./allK(logical(alls));
    pratio = vecmc(logical(allp))./allK(logical(allp));
    % non vector implementation for ADMAT
    for irxn = 1:length(Vind)     
        if size(sratio,2)>1
            sprod(irxn) = sratio(irxn,:).^-alls(logical(alls(:,irxn)),irxn);
            pprod(irxn) = pratio(irxn,:).^allp(logical(allp(:,irxn)),irxn);
        else
            sprod(irxn) = sratio(irxn).^-alls(logical(alls(:,irxn)),irxn);
            pprod(irxn) = pratio(irxn).^allp(logical(allp(:,irxn)),irxn);
        end
    end
    fwdflx = kfwd(Vind).*sprod;
    revflx = kbkw(Vind).*pprod;
%     fwdflx = kfwd(Vind).*prod(sratio.^-alls(logical(alls)),2);
%     revflx = kbkw(Vind).*prod(pratio.^allp(logical(allp)),2);    
    
    % set reverse flux for zero products = 0
    [~,rxn] = find(vecmc(logical(alls))==0);
    if ~isempty(rxn)
        revflx(rxn) = 0;
    end
    
    % set forwrd flux for zero substrate = 0
    [~,rxn] = find(vecmc(logical(allp))==0);
    if ~isempty(rxn)
        fwdflx(rxn) = 0;
    end
    
    % set reverse flux for irreversible reactions = 0
    revflx(~rev(Vind)) = 0;
    
    % numerator of kinetics
    nrflx = fwdflx-revflx;
    
    % denominator of kinetics
    % non vector implementation for ADMAT
    drflx = 1+sprod+pprod;
%     drflx = 1+...
%             prod(sratio.^-alls(logical(alls)),2)+...
%             prod(pratio.^allp(logical(allp)),2);
    
    % allosteric activation  
    % reaction that have an activation component
    [~,rxn] = find(SI(:,Vind)>0);
    Vact = Vind(rxn);
    allacs = SI(:,Vind);
    allacs(SI(:,Vind)<0) = 0;
    allKIa = KIact(:,Vind);
    acratio = vecmc(logical(allacs))./allKIa(logical(allacs));
    if ~isempty(acratio)
        idacflx = (0.5+(1-0.5).*acratio./(1+acratio)).^allacs(logical(allacs));
        for irxn = 1:size(acratio,1)
            if size(idacflx,2)>1
                acflx(Vact(irxn)) = prod(idacflx(irxn,:)');
            else
                acflx(Vact(irxn)) = idacflx(irxn);
            end
        end        
%         acflx(Vact) = prod((0.5+(1-0.5).*acratio./(1+acratio)).^allacs(logical(allacs)),2);
    else
        acflx(Vact) = 1;
    end
    nrflx(rxn) = acflx(Vact).*nrflx(rxn);
    
    % allosteric inhibition
    [~,rxn] = find(SI(:,Vind)<0);
    Vihb = Vind(rxn);
    allibs = SI(:,Vind);
    allibs(SI(:,Vind)>0) = 0;
    allKIi = KIihb(:,Vind);
    ibratio = vecmc(logical(allibs))./allKIi(logical(allibs));
    if ~isempty(ibratio)
        idibflx = (0.5+(1-0.5).*ibratio./(1+ibratio)).^allibs(logical(allibs));
        for irxn = 1:size(ibratio,1)
            if size(idibflx,2)>1
                ibflx(Vihb(irxn)) = prod(idibflx(irxn,:)');
            else
                ibflx(Vihb(irxn)) = idibflx(irxn);
            end
        end 
%         ibflx(Vihb) = prod((0.5+(1-0.5)*ibratio/(1+ibratio)).^-allibs(logical(allibs)),2);
    else
        ibflx(Vihb) = 1;
    end
    nrflx(rxn) = ibflx(Vihb).*nrflx(rxn);
    
    % other activation
%     spacratio = vecmc(SI(:,Vind)>0,Vind)./KIact(SI(:,Vind)>0,Vind);
%     if ~isempty(spacratio)
%         acdr(Vind) = sum((1/spacratio).^SI(SI(:,Vind)>0,Vind));
%     end
%     drflx = drflx+acdr(Vind);
%     
%     % other inhibition
%     spibratio = vecmc(SI(:,Vind)<0,Vind)./KIihb(SI(:,Vind)<0,Vind);
%     if ~isempty(spibratio)
%         ibdr(Vind) = sum(spibratio.^(-SI(SI(:,Vind)<0,Vind)));
%     end
%     drflx = drflx+ibdr(Vind);
    
    % scale flux
    vflux(Vind,ic) = scale_flux(nrflx./drflx);
    
    % scale all ---> A (no substrate) reactions to zero
%     [met,rxn] = find(S(S(:,Vind)<0,:)==0);
%     vflux(rxn,ic) = 0;
    
    % scale all A ---> (no product) reactions to zero 
%     [met,rxn] = find(S(S(:,Vind)>0,:)==0);
%     vflux(rxn,ic) = 0;
    
    % multiply with Vmax or enzyme concentration
    flux(Vind,ic) = Vmax(Vind,ic).*vflux(Vind,ic);    
end

flux = flux(Vind,:);
vflux = vflux(Vind,:);