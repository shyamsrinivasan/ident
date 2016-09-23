function [fluxbm] = BMKinetics(model,pvec,M,Vid)
[nvar,nc] = size(M);
S = model.S;
remid = model.remid;
nrxn = model.nt_rxn;
K = pvec.K;

fluxbm = zeros(nvar,1);


% Vid = bmrxn
% bmrxn = find(strcmpi(model.rxns,'biomass'));

ace = strcmpi(model.mets,'ac[e]');
ackr = strcmpi(model.rxns,'ACKr');
if nc>1
    alpha = M(ace,:)./repmat(K(ace,ackr),1,nc);
else
    alpha = M(ace,:)./K(ace,ackr);
end
alpha = alpha./(1+alpha);

for irxn = 1:nrxn
    if ~isempty(Vid)
        if ismember(irxn,Vid)    
            alls = S(:,irxn);allp = S(:,irxn);
            alls(S(:,irxn)>0) = 0;allp(S(:,irxn)<0) = 0;
            alls(remid,:) = 0;allp(remid,:) = 0;
            if nc>1
                fluxbm(logical(alls)) = -repmat(K(logical(alls),irxn),1,nc).*M(logical(alls),:);
            else
                fluxbm(logical(alls)) = -K(logical(alls),irxn).*M(logical(alls),:);
            end
            if nc>1
                fluxbm(logical(allp)) = repmat(K(logical(allp),irxn),1,nc).*M(logical(allp),:);
            else
                fluxbm(logical(allp)) = K(logical(allp),irxn).*M(logical(allp),:);
            end
            fluxbm = repmat(alpha,nvar,1).*fluxbm;
        end
    end
end
fluxbm = fluxbm./3600;



% he = find(strcmpi(model.mets,'h[e]'));
% hc = find(strcmpi(model.mets,'h[c]'));
% h2o = find(strcmpi(model.mets,'h2o[c]'));
% fdp = strcmpi(model.mets,'fdp[c]');
% pep = strcmpi(model.mets,'pep[c]');
% % vectorize mc
% vec_mc = repmat(mc,1,length(Vid));
% % vec_mc = sparse(vec_mc);
% vec_mc(~S(:,Vid)) = 0;
% 
% % vectorized mc for activation
% actvec_mc = repmat(mc,1,length(Vid));
% % actvec_mc = sparse(actvec_mc);
% actvec_mc(SI(:,Vid)~=1) = 0;
% 
% % vectorized mc for inhibition
% ihbvec_mc = repmat(mc,1,length(Vid));
% % ihbvec_mc = sparse(ihbvec_mc);
% ihbvec_mc(SI(:,Vid)~=-1) = 0;
% 
% % if rxn id is biomass reaction
% if Vid == model.bmrxn
%     ratio = 1+vec_mc(fdp)./K(fdp,Vid);
%     vflux(Vid,:) = (ratio-1).*(ratio).^3./(ratio.^4+4e6.*...
%                    (1+actvec_mc(pep)./KIact(pep,Vid)).^(-4));    
%     vflux(Vid) = scale_flux(vflux(Vid,:));
%     flux(Vid,:) = Vmax(Vid).*vflux(Vid,:);
% end 

% flux = flux(Vid,:);
% vflux = vflux(Vid,:);