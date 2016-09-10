function pvec = estimateKm(pvec,sbid,prid,mc,Ksbackup,Kpbackup,Vind,rerun)
% suprceded by samplesigma.m and samplekcat.m
if nargin<8
    rerun = 0;
end
sbfn = find(sbid);
prfn = find(prid);

sbzro = sbfn(mc(sbfn)==0);
przro = prfn(mc(prfn)==0);
    
nsb = length(find(sbid));
npr = length(find(prid));

Kscol = zeros(length(sbid),1);
Kpcol = zeros(length(prid),1);

for irxn = 1:length(Vind)
    %enzyme saturation
    sigma = random(makedist('Uniform'),...
                   nsb+npr,...
                   1);

    %determine kinetic parameters from concentration and sigma
    %substrates
    if ~isempty(find(sbid,1))
        sb_rat = sigma(1:nsb)./(1-sigma(1:nsb));
        Ksb = mc(logical(sbid))./sb_rat;
        if rerun
            Kscol(logical(sbid),1) = Ksbackup;
        else
            Kscol(logical(sbid),1) = pvec.K(logical(sbid),Vind(irxn));
        end
        if any(Kscol==1)
            pvec.K(Kscol==1,Vind(irxn))=Ksb(Kscol(logical(sbid),1)==1);
            % pvec.K(logical(sbid),Vind(irxn))
            pvec.K(sbzro,Vind(irxn)) = 1;
            pvec.Kind(Kscol==1,Vind(irxn))=1;
        end
    end

    %products
    if ~isempty(find(prid,1))
        pr_rat = sigma(nsb+1:nsb+npr)./(1-sigma(nsb+1:nsb+npr));
        Kpr = mc(logical(prid))./pr_rat;
        if rerun
            Kpcol(logical(prid),1) = Kpbackup;
        else
            Kpcol(logical(prid),1) = pvec.K(logical(prid),Vind(irxn));
        end
        if any(Kpcol==1)
            pvec.K(Kpcol==1,Vind(irxn))=Kpr(Kpcol(logical(prid),1)==1);
            pvec.K(przro,Vind(irxn)) = 1;
            pvec.Kind(Kpcol==1,Vind(irxn))=1;
        end
    end
end