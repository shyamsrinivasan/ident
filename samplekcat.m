function kcatbkw = samplekcat(delGlb,delGub,delG,rdir,sbid,prid,mc)

%sample delGr
%positive direction
if rdir > 0
    delG = random(makedist('Uniform','lower',delGlb,'upper',0),1,1);
elseif rdir < 0
    delG = random(makedist('Uniform','lower',0,'upper',delGub),1,1);
end

%reverse kcat for reversible reactions
%normalized substrates
nrm_sb = prod(mc(sbid)./K(sbid,irxn));

%normalized product
nrm_pr = prod(mc(prid)./K(prid,irxn));

kcatbkw = kcatfwd*nrm_sb/nrm_pr*exp(delG/RT);
    
    

