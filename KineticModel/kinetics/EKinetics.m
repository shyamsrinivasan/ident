function [flux,vflux] = EKinetics(model,pvec,M,VFex)
[~,nc] = size(M);
S = model.S;
nrxn = model.nt_rxn;
% rev = model.rev;
remid = model.remid;
D = model.D;

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
            fwdflx = D.*thetas; % mmole/Lc/h
            revflx = zeros(1,nc);
            if any(Vuptake(irxn))
                % assume feeds of Vuptake mmole/Lc
                revflx = D.*0.5; % repmat(Vuptake(irxn),1,nc); % mmole/Lc/h
%                 revflx = repmat(Vuptake(irxn),1,nc);
            end        
            nrflx = fwdflx-revflx;
            drflx = ones(1,nc);
            vflux(irxn,:) = scale_flux(nrflx./drflx);
            flux(irxn,:) = Vmax(irxn).*vflux(irxn,:);
            % convert mmole/Lc/h -> mmole/Lc/s
            flux(irxn,:) = flux(irxn,:)./repmat(3600,1,nc);
        end
    end
end
flux = flux(VFex,:);
vflux = vflux(VFex,:);
