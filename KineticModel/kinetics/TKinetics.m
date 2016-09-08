function [flux,vflux] = TKinetics(model,pvec,M,Vex)
[~,nc] = size(M);
S = model.S;
nrxn = model.nt_rxn;
rev = model.rev;
K = pvec.K;
kfwd = pvec.kfwd;
kbkw = pvec.krev;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,nc);
flux = zeros(nrxn,nc);

% vflux = cons(vflux,M);
% flux = cons(flux,M);

% vecmc = repmat(M,1,nrxn);

pie = strcmpi(model.mets,'pi[e]');
% eliminate consideration for excess cofators
% pi[c],pi[e],h[c],h[e],h2o[c]
he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
S([he hc h2o],:) = 0;

% h2o = find(strcmpi(model.mets,'h2o[c]'));

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

%     %kinetics - substrate and product
%     sbid = logical(model.S(:,Vex(irxn))<0);
%     prid = logical(model.S(:,Vex(irxn))>0); 
%     
%     % remove protons
%     sbid([he hc]) = 0;
%     prid([he hc]) = 0;
%     
%     if any(strcmpi(model.rxns{Vex(irxn)},'O2t'))
%         nr_flux = kfwd(Vex(irxn))*prod(M(sbid)./K(sbid,Vex(irxn)))/...
%                   (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
%                   prod(M(prid)./K(prid,Vex(irxn))));
%     elseif any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%         Kapie = 0.89; % mM
%         nr_flux = kfwd(Vex(irxn))*prod(M(sbid)./K(sbid,Vex(irxn)))/...
%                   (1+prod(M(sbid)./K(sbid,Vex(irxn)))+...
%                   prod(M(prid)./K(prid,Vex(irxn))))*...
%                   1/(1+Kapie/M(sbid));

    
%     if any(strcmpi(model.rxns{Vex(irxn)},'o2t'))
% %         Do2 = 2.1e-9; %m2/s
% %         Acell = 4.42e-12; %m2
% %         vflux(Vex(irxn)) = Do2/Acell*(mc(sbid)-mc(prid));
%         vflux(Vex(irxn)) = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
%                            (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
%                            prod(mc(prid)./K(prid,Vex(irxn))));
% %         Vmax(Vex(irxn)) = 1;
%     end

for irxn = 1:nrxn
    if ismember(irxn,Vex)
        alls = S(:,irxn);allp = S(:,irxn);
        alls(S(:,irxn)>0) = 0;allp(S(:,irxn)<0) = 0;
        sratio = M(logical(alls),:)./K(logical(alls),irxn);
        pratio = M(logical(allp),:)./K(logical(allp),irxn);
        thetas = prod(sratio.^-alls(logical(alls)),1);
        thetap = prod(pratio.^allp(logical(allp)),1);
        fwdflx = kfwd(irxn).*thetas;
        revflx = kbkw(irxn).*thetap;
        % set reverse flux for irreversible reactions = 0
        if ~rev(irxn)
            revflx = zeros(1,length(revflx));
        end
        % numerator of kinetics
        nrflx = fwdflx-revflx;
        % denominator of kinetics
        drflx = 1+thetas+thetap;
        % scale flux
        vflux(irxn,:) = scale_flux(nrflx./drflx);
        
        % fluxes for O2t and PIt2r from Vex - overwrite nrflx and drflx
        if irxn == find(strcmpi(model.rxns,'O2t'))
        end
        if irxn == find(strcmpi(model.rxns,'PIt2r'))
            Kapie = 0.89; % mM
            vflux(irxn,:) =...
            scale_flux(fwdflx./drflx.*1./(1+Kapie./M(pie,:)));
        end
        flux(irxn,:) = Vmax(irxn).*vflux(irxn,:); 
    end
end
flux = flux(Vex,:);
vflux = vflux(Vex,:);