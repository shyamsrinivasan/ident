function [newkf,newkr] = samplekcat(alls,allp,rev,kfwd,krev,K,Keq)
nmodels = size(K,2);
% pvec = samplekcat(model,pvec,sbid,prid,irxn,mc,kfwdbkup,kbkwbkup,rerun)
% [newkf,newkr] = samplekcat(model,Vind,kfwd,krev,K,Keq,nmodels)

% if nargin<6
%     nmodels = 1;
% end

% nrxn = model.nt_rxn;
% S = model.S;

% R = 0.008314; %kJ/mol.K
% T = 298; %K
% RT = R*T;

sstoich = repmat(-alls(logical(alls)),1,nmodels);    
pstoich = repmat(allp(logical(allp)),1,nmodels);
Keq = repmat(Keq,1,nmodels);

if ~rev
    if ~isnan(kfwd)
        newkf = repmat(kfwd,1,nmodels);
    end
    newkr = zeros(1,nmodels);
elseif rev
    if ~isnan(kfwd)
        newkf = repmat(kfwd,1,nmodels);
    end
    if ~isnan(krev)
        newkr = repmat(krev,1,nmodels);
    end
end

if isnan(kfwd)
    krbk = repmat(krev,1,nmodels);    
    % new kfwd as nmodels vector
    newkf = Keq.*krbk.*prod(K(1:length(find(alls)),:).^sstoich,1)./...
            prod(K(length(find(alls))+1:length(find(alls))+length(find(allp)),:).^pstoich,1);
    newkr = krbk;
end

if isnan(krev)
    kfbk = repmat(kfwd,1,nmodels);
    % new krev as nmodels vector
    newkr = kfbk./Keq.*...
    prod(K(length(find(alls))+1:length(find(alls))+length(find(allp)),:).^pstoich,1)./...
    prod(K(1:length(find(alls)),:).^sstoich,1);
    newkf = kfbk;
end



%     pvec.kfwd(irxn) = model.Keq(irxn)*kcatbkw*...
%                           prod(K(sbid,irxn).^Sb)/prod(K(prid,irxn).^Sp);
%     pvec.kfwd(irxn) = nrm_pr/nrm_sb*kcatbkw*exp(-delGr/RT);    

%     pvec.krev(irxn) = (kcatfwd/model.Keq(irxn))*...
%                           prod(K(prid,irxn).^Sp)/prod(K(sbid,irxn).^Sb);
%     pvec.krev(irxn) = nrm_sb/nrm_pr*kcatfwd*exp(delGr/RT);

% if model.Vss(irxn)==0
%     pvec.kfwd(irxn) = 0;
%     pvec.krev(irxn) = 0;
% end
% fprintf('kcat+ = %3.6g \t kcat- = %3.6g \t Keq = %3.6g\t',pvec.kfwd(irxn),...
%                                          pvec.krev(irxn),...
%                                          model.Keq(irxn));
% fprintf('vflux = %3.6g\n',...                                   
% pvec.kfwd(irxn)*nrm_sb - pvec.krev(irxn)*nrm_pr);



    
    

