function pvec = samplekcatTKinetic(model,pvec,mc)

R = 0.008314; %kJ/mol.K
T = 298; %K
RT = R*T;

S = model.S;

h2o = find(strcmpi(model.mets,'h2o[c]'));
pic = find(strcmpi(model.mets,'pi[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));

for irxn = 1:length(Vex)
    sbid = S(:,Vind(irxn))<0;    
    prid = S(:,Vind(irxn))>0;
    
    Kscol = zeros(length(find(sbid)),1);
    Kpcol = zeros(length(find(prid)),1);
    
    if ~strcmpi(model.rxns{Vex(irxn)},'h2ot')
        sbid(h2o) = 0;
        prid(h2o) = 0;
    end
    if ~strcmpi(model.rxns{Vex(irxn)},'PIt2r')
        sbid([pic pie hc he]) = 0;
        prid([pic pie hc he]) = 0;
    else
        sbid([hc he]) = 0;
        prid([hc he]) = 0;
    end
    
    sbfn = find(sbid);
    prfn = find(prid);
    
    sbzro = sbfn(mc(sbfn)==0);
    przro = prfn(mc(prfn)==0);
    
    nsb = length(find(sbid));
    npr = length(find(prid));
    
    %enzyme saturation
    sigma = random(makedist('Uniform'),...
                   nsb+npr,...
                   1);
               
    if ~isempty(find(sbid,1))
        sb_rat = sigma(1:nsb)./(1-sigma(1:nsb));
        Ksb = mc(logical(sbid))./sb_rat;
        Kscol(logical(sbid),1) = pvec.K(logical(sbid),Vind(irxn));
        if any(Kscol==1)
            pvec.K(Kscol==1,Vind(irxn))=Ksb(pvec.K(logical(sbid),Vind(irxn))==1);
            pvec.K(sbzro,Vind(irxn)) = 1;
            pvec.Kind(Kscol==1,Vind(irxn))=1;
        end
    end
    
    if ~isempty(find(prid,1))
        pr_rat = sigma(nsb+1:nsb+npr)./(1-sigma(nsb+1:nsb+npr));
        Kpr = mc(logical(prid))./pr_rat;
        Kpcol(logical(prid),1) = pvec.K(logical(prid),Vind(irxn));
        if any(Kpcol==1)
            pvec.K(Kpcol==1,Vind(irxn))=Kpr(pvec.K(logical(prid),Vind(irxn))==1);
            pvec.K(przro,Vind(irxn)) = 1;
            pvec.Kind(Kpcol==1,Vind(irxn))=1;
        end
    end
    
    if isnan(pvec.kcat_fwd(Vex(irxn)))
        pvec.kcat_fwd(Vex(irxn)) = pvec.kcat_bkw(Vex(irxn))*1*...
                                   prod(pvec.K(sbid,Vex(irxn)))/...
                                   prod(pvec.K(prid,Vex(irxn)));
    end
    if isnan()
        pvec.kcat_bkw(Vex(irxn)) = (pvec.kcat_fwd(Vex(irxn))/1)*...
                                    prod(pvec.K(prid,Vex(irxn)))/...
                                    prod(pvec.K(sbid,Vex(irxn)));
    end
        
        
    
    
end