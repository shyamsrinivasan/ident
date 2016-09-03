function [flux,vflux] = TKinetics(model,pvec,M,Vex)
Vmax = pvec.Vmax;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
K = pvec.K;
S = model.S;
rev = model.rev;

nrxn = model.nt_rxn;
flux = zeros(nrxn,1);
flux = cons(flux,M);
vflux = zeros(nrxn,1);
vflux = cons(vflux,M);

he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
% h2o = find(strcmpi(model.mets,'h2o[c]'));
pie = strcmpi(model.mets,'pi[e]');
% pic = find(strcmpi(model.mets,'pi[c]'));
% co2 = find(strcmpi(model.mets,'co2[c]'));
    
% q8 = find(strcmpi(model.mets,'q8[c]'));
% q8h2 = find(strcmpi(model.mets,'q8h2[c]'));

% vmet = [he hc pie pic h2o co2];   

% vpts = find(strcmpi(model.rxns,'GLCpts'));
% vthd2 = find(strcmpi(model.rxns,'THD2'));%nadph --->
% vatps = find(strcmpi(model.rxns,'ATPS4r'));%atp ---> adp
% vnad16 = find(strcmpi(model.rxns,'NADH16'));%nad ---> nadh
% vcyt = find(strcmpi(model.rxns,'CYTBD'));%o2 --->
% vspl = [vpts vthd2 vatps vnad16 vcyt];

% nadph = strcmpi(model.mets,'nadph[c]');
% nadp = strcmpi(model.mets,'nadp[c]');
% nad = strcmpi(model.mets,'nad[c]');
% nadh = strcmpi(model.mets,'nadh[c]');
% atp = strcmpi(model.mets,'atp[c]');
% adp = strcmpi(model.mets,'adp[c]');
% o2 = strcmpi(model.mets,'o2[c]');

for irxn = 1:length(Vex)    
        
    %kcat
%     kfwd = pvec.kcat_fwd(Vex(irxn));
%     kbkw = pvec.kcat_bkw(Vex(irxn));

    nmet = size(S,1);
    %kinetics - substrate and product
    sbid = logical(model.S(:,Vex(irxn))<0);
    prid = logical(model.S(:,Vex(irxn))>0); 
    
    % remove protons
    sbid([he hc]) = 0;
    prid([he hc]) = 0;
    
    if any(strcmpi(model.rxns{Vex(irxn)},'O2t'))
        nr_flux = kfwd(Vex(irxn))*prod(M(sbid)./K(sbid,Vex(irxn)))/...
                  (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
                  prod(M(prid)./K(prid,Vex(irxn))));
    elseif any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
        Kapie = 0.89; % mM
        nr_flux = kfwd(Vex(irxn))*prod(M(sbid)./K(sbid,Vex(irxn)))/...
                  (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
                  prod(M(prid)./K(prid,Vex(irxn))))*...
                  1/(1+Kapie/M(sbid));
    else
        if rev(Vex(irxn))
            nr_flux = kfwd(Vex(irxn))*(prod(M(sbid)./K(sbid,Vex(irxn)))-...
                      prod(M(prid)./K(prid,Vex(irxn))))/...
                      (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
                      prod(M(prid)./K(prid,Vex(irxn))));
        elseif ~rev(Vex(irxn))
            nr_flux = kfwd(Vex(irxn))*prod(M(sbid)./K(sbid,Vex(irxn)))/...
                      (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
                      prod(M(prid)./K(prid,Vex(irxn))));
        end
    end
    if any(sbid) && any(prid)
        vflux(Vex(irxn)) = scale_flux(nr_flux);
    else
        vflux(Vex(irxn)) = 0;
    end
 
%     if model.rev(Vex(irxn))         
%         if all(mc(sbid)>0) && all(mc(prid)>0)
%             nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn))) -...
%                       mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
%         elseif all(mc(sbid)>0)
%             nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));
%         elseif all(mc(prid)>0)
%             nr_flux = -mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
%         end
%         if any(sbid) && any(prid)
%             % Denominator - 1.6
%             dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));            
%             % dr_pr
%             dr_pr = 1+mc(prid)./K(prid,Vex(irxn));            
%         else
%             dr_sb = 0;
%             dr_pr = 0;
%         end        
%     elseif ~model.rev(Vex(irxn))        
%         nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));  
%         if any(sbid) && any(prid)
%             dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));            
%         else
%             dr_sb = 0;
%         end
%     end
%     
%     if any(sbid) && any(prid)
%         %Denominator - 1.6
%         dr_flux = prod(dr_sb)+prod(dr_pr)-1;
%         vflux(Vex(irxn)) = scale_flux(nr_flux/dr_flux);
%     else
%         vflux(Vex(irxn)) = 0;
%     end
    
%     if any(strcmpi(model.rxns{Vex(irxn)},'o2t'))
% %         Do2 = 2.1e-9; %m2/s
% %         Acell = 4.42e-12; %m2
% %         vflux(Vex(irxn)) = Do2/Acell*(mc(sbid)-mc(prid));
%         vflux(Vex(irxn)) = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
%                            (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
%                            prod(mc(prid)./K(prid,Vex(irxn))));
% %         Vmax(Vex(irxn)) = 1;
%     end
%     flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
    vflux(Vex(irxn)) = scale_flux(vflux(Vex(irxn)));
    flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
end

% partial vector implementation to reduce overhead and use with ADMAT
sprod = ones(length(Vex),1);
pprod = ones(length(Vex),1);
vecmc = repmat(M,1,length(Vex));
alls = S(:,Vex);allp = S(:,Vex);
alls(S(:,Vex)>0) = 0;allp(S(:,Vex)<0) = 0;
allK = K(:,Vex);

sratio = vecmc(logical(alls))./allK(logical(alls));
pratio = vecmc(logical(allp))./allK(logical(allp));
for irxn = 1:length(Vex)
    if size(sratio,2)>1
        sprod(irxn) = prod(sratio(irxn,:));
        pprod(irxn) = prod(pratio(irxn,:));
    else
        sprod(irxn) = sratio(irxn);
        pprod(irxn) = pratio(irxn);
    end
end
fwdflx = kfwd(Vex).*sprod;
revflx = kfwd(Vex).*pprod;

% set reverse flux for zero products = 0
[~,rxn] = find(vecmc(logical(allp))==0);
if ~isempty(rxn)
    revflx(rxn) = 0;
end

% set forwrd flux for zero substrate = 0
[~,rxn] = find(vecmc(logical(alls))==0);
if ~isempty(rxn)
    fwdflx(rxn) = 0;
end

% set reverse flux for irreversible reactions = 0
revflx(~rev(Vex)) = 0;

% numerator of kinetics
nrflx = fwdflx-revflx;

% denominator of kinetics
% non vector implementation for ADMAT
drflx = 1+sprod+pprod;

% scale flux
vflux(Vex) = scale_flux(nrflx./drflx);

% fluxes for O2t and PIt2r from Vex - overwrite nrflx and drflx
tfo2 = ismember(find(strcmpi(model.rxns,'O2t')),Vex); 
if ~isempty(tfo2) && any(tfo2)
    vflux(tfo2) = scale_flux(fwdflx(tfo2)/drflx(tfo2));
end

tfpit2r = ismember(find(strcmpi(model.rxns,'PIt2r')),Vex); 
if ~isempty(tfpit2r) && any(tfpit2r)
    Kapie = 0.89; % mM
    vflux(tfpit2r) =...
    scale_flux(fwdflx(tfpit2r)./drflx(tfpit2r).*1./(1+Kapie./vecmc(pie,tfpit2r)));
end
        
flux(Vex) = Vmax(Vex).*vflux(Vex);  

flux = flux(Vex);
vflux = vflux(Vex);
