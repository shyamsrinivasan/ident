function flux = EKinetics(model,pvec,M,VFex)
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
        sratio = M(logical(alls),:);
        thetas = sratio;
        fwdflx = 0.1*thetas;
        revflx = 0;
        if any(Vuptake(irxn))
            revflx = Vuptake(irxn);
        end        
        nrflx = fwdflx-revflx;
        drflx = 1;
        vflux(irxn,:) = scale_flux(nrflx./drflx);
        flux(irxn,:) = Vmax(irxn).*vflux(irxn,:);
    end
end
flux = flux(VFex);
vflux = vflux(VFex);
