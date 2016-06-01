function pvec = samplekcat(model,pvec,sbid,prid,irxn,mc,kfwdbkup,kbkwbkup,rerun)

% delGr = pvec.delGr(irxn);     
K = pvec.K;

R = 0.008314; %kJ/mol.K
T = 298; %K
% RT = R*T;

%normalized substrates
% nrm_sb = prod(mc(sbid)./K(sbid,irxn));

%normalized product
% nrm_pr = prod(mc(prid)./K(prid,irxn));

Sb = -model.S(sbid,irxn); 
Sp = model.S(prid,irxn);   

if rerun
    kcatfwd = kfwdbkup;
    kcatbkw = kbkwbkup;
else
    kcatfwd = pvec.kcat_fwd(irxn);
    kcatbkw = pvec.kcat_bkw(irxn);
end

if isnan(kcatfwd)
    pvec.kcat_fwd(irxn) = model.Keq(irxn)*kcatbkw*...
                          prod(K(sbid,irxn).^Sb)/prod(K(prid,irxn).^Sp);
%     pvec.kcat_fwd(irxn) = nrm_pr/nrm_sb*kcatbkw*exp(-delGr/RT);    
end
if isnan(kcatbkw)
    pvec.kcat_bkw(irxn) = (kcatfwd/model.Keq(irxn))*...
                          prod(K(prid,irxn).^Sp)/prod(K(sbid,irxn).^Sb);
%     pvec.kcat_bkw(irxn) = nrm_sb/nrm_pr*kcatfwd*exp(delGr/RT);
end
if model.Vss(irxn)==0
    pvec.kcat_fwd(irxn) = 0;
    pvec.kcat_bkw(irxn) = 0;
end
% fprintf('kcat+ = %3.6g \t kcat- = %3.6g \t Keq = %3.6g\t',pvec.kcat_fwd(irxn),...
%                                          pvec.kcat_bkw(irxn),...
%                                          model.Keq(irxn));
% fprintf('vflux = %3.6g\n',...                                   
% pvec.kcat_fwd(irxn)*nrm_sb - pvec.kcat_bkw(irxn)*nrm_pr);



    
    

