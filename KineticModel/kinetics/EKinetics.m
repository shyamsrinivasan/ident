function [flux,vflux] = EKinetics(model,pvec,M,VFex)
[~,nc] = size(M);
S = model.S;
nrxn = model.nt_rxn;
% rev = model.rev;
remid = model.remid;

% K = pvec.K;
% kfwd = pvec.kfwd;
% kbkw = pvec.krev;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,nc);
flux = zeros(nrxn,nc);
% DVX = zeros(length(M),nrxn);

if isfield(model,'Vuptake')
    Vuptake = model.Vuptake;
end

for irxn = 1:nrxn
    if ismember(irxn,VFex)
        alls = S(:,irxn);alls(S(:,irxn)>0) = 0;alls(remid,:) = 0;
        if any(alls)
            sratio = M(logical(alls),:);
            thetas = sratio;
            fwdflx = 0.1.*thetas;
            revflx = zeros(1,nc);
            if any(Vuptake(irxn))
                revflx = repmat(Vuptake(irxn),1,nc);
            end        
            nrflx = fwdflx-revflx;
            drflx = ones(1,nc);
            vflux(irxn,:) = scale_flux(nrflx./drflx);
            flux(irxn,:) = Vmax(irxn).*vflux(irxn,:)./repmat(3600,1,nc);
        end
    end
end
flux = flux(VFex,:);
vflux = vflux(VFex,:);
