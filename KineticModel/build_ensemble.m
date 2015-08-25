%Sample Quantities to build ensemble of models
function [ensb,kflag,Vflag] = build_ensemble(nmodels,model,pvec,MC,KVl)

nsamples = nmodels;

%glcpts
Vind = [model.Vind find(strcmpi(model.rxns,'GLCpts'))];
Vex = setdiff(model.Vex,Vind);

nrxn = length(Vind);
nt_rxn = model.nt_rxn;
ntmetab = size(model.S,1);
for isample = 1:nsamples
    mname = sprintf('model%d',isample);
    ensb.(mname).K = pvec.K;
    ensb.(mname).Kind = sparse(ntmetab,nt_rxn);
    ensb.(mname).KIact = pvec.KIact;
    ensb.(mname).KIihb = pvec.KIihb;
    ensb.(mname).KAind = sparse(ntmetab,nt_rxn);
    ensb.(mname).KIind = sparse(ntmetab,nt_rxn);
    ensb.(mname).Vmax = pvec.Vmax;
    ensb.(mname).KVl = KVl;
    ensb.(mname).kcat_ratio = zeros(nt_rxn,1);
    ensb.(mname).kcat_fwd = zeros(nt_rxn,1); 
    ensb.(mname).kcat_bkw = zeros(nt_rxn,1);
end   

  %Remove kinetic consideration for mets like water which are
  %non-limiting always
for isample = 1:nsamples
    mname = sprintf('model%d',isample);  
    jmetab = 0;
    kreg = 0;   
    Vmax = zeros(nt_rxn,1);
    for irxn = 1:nrxn  
        sub = find(model.S(:,Vind(irxn))<0);
        prd = find(model.S(:,Vind(irxn))>0);
        nmet = length(sub)+length(prd);       
        Kscol = zeros(ntmetab,1);
        Kpcol = zeros(ntmetab,1);
        nmetab = length(sub);
                
        %Substrate/Product Samples   
        met_sat = rand(nmet,1);%uniform pseudo random saturations (0-1)
%         met_sat = betarnd(1.5,4.5,nmetab(irxn),1); %Beta Distribution
        nsubs = length(sub);        
        if nsubs ~= 0
            sub_ratio = met_sat(1:nsubs)./(1-met_sat(1:nsubs));
            Ksubs = MC(sub)./sub_ratio;
            Kscol(sub,1) = ensb.(mname).K(sub,Vind(irxn));
            if any(Kscol==1)
                ensb.(mname).K(Kscol==1,Vind(irxn)) =...
                Ksubs(ensb.(mname).K(sub,Vind(irxn))==1); 
                ensb.(mname).Kind(Kscol==1,Vind(irxn)) = 1;
            end
        end
        if nmet-nsubs ~= 0
            prod_ratio = met_sat(nsubs+1:nmet)./(1-met_sat(nsubs+1:nmet));
            Kprod = MC(prd)./prod_ratio;
            Kpcol(prd,1) = ensb.(mname).K(prd,Vind(irxn));
            if any(Kpcol == 1)
                ensb.(mname).K(Kpcol==1,Vind(irxn)) =...
                Kprod(ensb.(mname).K(prd,Vind(irxn))==1);
                ensb.(mname).Kind(Kpcol==1,Vind(irxn)) = 1;
            end
        end
        
        act = find(model.SI(:,Vind(irxn))>0);
        ihb = find(model.SI(:,Vind(irxn))<0);
        nreg = length(act)+length(ihb);
        KAcol = zeros(ntmetab,1);
        KIcol = zeros(ntmetab,1);   
        
        %Regulator Samples 
        reg_sat = rand(nreg,1);%uniform pseudo random saturations (0-1)
%         reg_sat = betarnd(1.5,4.5,nreg(irxn),1);
        nact = length(act);        
        %Activation
        if nact ~= 0            
            act_ratio = reg_sat(1:nact)./(1-reg_sat(1:nact));
            KAreg = MC(act)./act_ratio;
            KAcol(act,1) = ensb.(mname).KIact(act,Vind(irxn));
            if any(KAcol == 1)
                ensb.(mname).KIact(KAcol==1,Vind(irxn))=...
                KAreg(ensb.(mname).KIact(act,Vind(irxn))==1);
                ensb.(mname).KAind(KAcol==1,Vind(irxn)) = 1;
            end
        end
        
        %Inhibition
        if nreg-nact ~= 0
            ihb_ratio = reg_sat(nact+1:nreg)./(1-reg_sat(nact+1:nreg));
            KIreg = MC(ihb)./ihb_ratio;
            KIcol(ihb,1) = ensb.(mname).KIihb(ihb,Vind(irxn));
            if any(KIcol==1)
                ensb.(mname).KIihb(KIcol==1,Vind(irxn))=...
                KIreg(ensb.(mname).KIihb(ihb,Vind(irxn))==1);
                ensb.(mname).KIind(KIcol==1,Vind(irxn)) = 1;
            end
        end
        jmetab = jmetab + nmetab;
        kreg = kreg + nreg;
%         ***********================**************
%         if ~isnan(model.Kcat(Vind(irxn)))
%             [kcat_ratio,kcat_fwd] = haldane_kcat(model,model.Kcat,...
%                                                  sub,prd,...
%                                                  Vind(irxn),ensb.(mname));
%         else%use sampled KVl
            [kcat_ratio,kcat_fwd,kcat_bkw] = haldane_kcat(model,KVl,...
                                                 sub,prd,...
                                                 Vind(irxn),ensb.(mname));
            
%         end
        ensb.(mname).kcat_ratio(Vind(irxn)) = kcat_ratio;
        ensb.(mname).kcat_fwd(Vind(irxn)) = kcat_fwd;   
        ensb.(mname).kcat_bkw(Vind(irxn)) = kcat_bkw;
        [~,vflux] = ConvinienceKinetics(model,ensb.(mname),MC,Vind(irxn));
        if model.Vss(Vind(irxn)) ~= 0
            Vmax(Vind(irxn)) = model.Vss(Vind(irxn))/vflux(Vind(irxn));
        else
            Vmax(Vind(irxn)) = 0;
        end
    end   
%     [~,vflux] = ConvinienceKinetics(model,ensb.(mname),MC,model.bmrxn,Vind);
    oldVmax = ensb.(mname).Vmax;
    ensb.(mname).Vmax(Vind) = Vmax(Vind);%model.Vss(Vind)./vflux(Vind);
    ensb.(mname).Vmax(setdiff(1:nt_rxn,Vind)) = 1;
    if any(~isnan(oldVmax))
        ensb.(mname).Vmax(~isnan(oldVmax)) = oldVmax(~isnan(oldVmax));
    end    
    %exhcnage reaction kcat_fwd and kcat_bkw
    for irxn = 1:length(Vex)
        if ~isnan(model.Kcat(Vex(irxn)))
            ensb.(mname).kcat_fwd(Vex(irxn)) = model.Kcat(Vex(irxn));
        %kla for o2 transfer - OTR
            if Vex(irxn) == find(strcmpi(model.rxns,'O2t'))
                ensb.(mname).kcat_bkw(Vex(irxn)) = model.Kcat(Vex(irxn));
            else
                ensb.(mname).kcat_bkw(Vex(irxn)) = model.Kcat(Vex(irxn));
            end
        else
           ensb.(mname).kcat_fwd(Vex(irxn)) = 0;
        end           
    end
end

%Select models from ensemble
kflag = 0;
Vflag = 0;
for is = 1:nsamples
    mname = sprintf('model%d',is);  
    kcatfl = zeros(length(Vind),1);
    for irxn = 1:length(Vind)
        if model.Vss(Vind(irxn)) < 0
            if ensb.(mname).kcat_bkw(Vind(irxn)) > ensb.(mname).kcat_fwd(Vind(irxn))
                kcatfl(irxn) = 1;
            end
        elseif model.Vss(Vind(irxn))>0
            if ensb.(mname).kcat_bkw(Vind(irxn)) < ensb.(mname).kcat_fwd(Vind(irxn))
                kcatfl(irxn) = 1;
            end
        else
            if ensb.(mname).kcat_fwd(Vind(irxn)) == 0
                kcatfl(irxn) = 1;
            end
        end
    end
    if all(kcatfl)
        kflag = 1;
    end
    if any(ensb.(mname).Vmax < 0)
        Vflag = 1;        
    end
end

% for is = 1:nsamples
%     mname = sprintf('model%d',is);  
%     if any(ensb.(mname).Vmax < 0)
%         flag = 1;        
%     end
% end