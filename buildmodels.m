function pvec = buildmodels(model,pvec,mc)

%reactions to consider for kinetics other than Vind
Vind = [model.Vind...
        find(strcmpi(model.rxns,'GLCpts'))];
%         find(strcmpi(model.rxns,'THD2'))...
%         find(strcmpi(model.rxns,'NADH16'))...
%         find(strcmpi(model.rxns,'ATPS4r'))];

Vind = setdiff(Vind,find(strcmpi(model.rxns,'ATPM')));

%metabolites that do not affect thermodynamic equilibrium   
vmet = [find(strcmpi(model.mets,'h[c]'))...
        find(strcmpi(model.mets,'h[e]'))...
        find(strcmpi(model.mets,'pi[c]'))...
        find(strcmpi(model.mets,'pi[e]'))...
        find(strcmpi(model.mets,'h2o[c]'))];
% q8 = find(strcmpi(model.mets,'q8[c]'));
% q8h2 = find(strcmpi(model.mets,'q8h2[c]'));

nrxn = length(Vind);
ntrxn = model.nt_rxn;
ntmet = model.nt_metab;

%backup known parameters from pvec
newp = struct();
newp.K = pvec.K;
newp.Kind = sparse(ntmet,ntrxn);
newp.Vmax = pvec.Vmax;
newp.kfwd = pvec.kcat_fwd;
newp.kbkw = pvec.kcat_bkw;

pvec.Kin = sparse(ntmet,ntrxn);
S = model.S;
check = zeros(ntrxn,1);
flux = zeros(ntrxn,1);

for irxn = 1:nrxn
    sbid = S(:,Vind(irxn))<0;    
    prid = S(:,Vind(irxn))>0;
    
    Kscol = zeros(length(find(sbid)),1);
    Kpcol = zeros(length(find(prid)),1);
    
    %no parameters for cofactors - assumed abundant 
    sbid(vmet) = 0;
    prid(vmet) = 0;
    
    nsb = length(find(sbid));
    npr = length(find(prid));
    
    %enzyme saturation
    sigma = random(makedist('Uniform'),...
                   nsb+npr,...
                   1);
    
    %determine kinetic parameters from concentration and sigma
    %substrates
    if ~isempty(find(sbid,1))
        sb_rat = sigma(1:nsb)./(1-sigma(1:nsb));
        Ksb = mc(logical(sbid))./sb_rat;
        Kscol(logical(sbid),1) = pvec.K(logical(sbid),Vind(irxn));
        if any(Kscol==1)
            pvec.K(Kscol==1,Vind(irxn))=Ksb(pvec.K(logical(sbid),Vind(irxn))==1);
            pvec.Kind(Kscol==1,Vind(irxn))=1;
        end
    end
    
    %products
    if ~isempty(find(prid,1))
        pr_rat = sigma(nsb+1:nsb+npr)./(1-sigma(nsb+1:nsb+npr));
        Kpr = mc(logical(prid))./pr_rat;
        Kpcol(logical(prid),1) = pvec.K(logical(prid),Vind(irxn));
        if any(Kpcol==1)
            pvec.K(Kpcol==1,Vind(irxn))=Kpr(pvec.K(logical(prid),Vind(irxn))==1);
            pvec.Kind(Kpcol==1,Vind(irxn))=1;
        end
    end
    
    %forward and backward catalytic rates
    %kfwd and kbkw
    %kfwd or kbkw is sampled basedon reaction directionality from FBA for
    %thermodynamic consistency
    %sampling done only for unknown values
    fprintf('irxn = %d \t delG = %3.6g\n',Vind(irxn),...
             pvec.delGr(Vind(irxn)));    
    pvec = samplekcat(pvec,sbid,prid,Vind(irxn),mc);
    
    pvec.Vmax(Vind(irxn)) = 1;
    
    %#check for vss and delGr direction    
    flux(Vind(irxn)) = CKinetics(model,pvec,mc,Vind(irxn));
    if pvec.delGr(Vind(irxn)) ~= 0
        if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<0
            check(Vind(irxn)) = 1;
        else
            check(Vind(irxn)) = -1;
        end    
    else
        if flux(Vind(irxn))*pvec.delGr(Vind(irxn))<1e-6
            check(Vind(irxn)) = 1;
        else
            check(Vind(irxn)) = -1;
        end
    end    
end

%restore Vmax from backup
pvec.Vmax = newp.Vmax;

%estimate Vmax
if all(check(Vind)>0)
    pvec.Vmax(pvec.delGr==0) = 0;
    pvec = findVmax(model,pvec,mc);
    
    %simple vmax = vss/ck
    for irxn = 1:length(Vind)
        [~,ck] = CKinetics(model,pvec,mc,Vind(irxn));
        if ck
            pvec.Vmax(Vind(irxn)) = model.Vss(Vind(irxn))/ck;
        else
            pvec.Vmax(Vind(irxn)) = 1;
        end
    end
else
    fprintf('Thermodynamically infeasible parameters\n');
    fprintf('Discontinuing\n');
    return
end
pvec.kcat_fwd(model.VFex) = 0;
pvec.kcat_bkw(model.VFex) = 0;

%check - calculate initial flux
flux = iflux(model,pvec,mc,flux);

    
    
            



        