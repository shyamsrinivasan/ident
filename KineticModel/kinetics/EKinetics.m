function flux = EKinetics(model,pvec,M,VFex)
[~,nc] = size(M);
flux = zeros(model.nt_rxn,nc);
% flux = cons(flux,M);

kfwd = pvec.kfwd;
% kcatbkw = pvec.kcat_bkw;
% flux = zeros(model.nt_rxn,1);
% model.Vuptake = zeros(model.nt_rxn,1);
% model.Vuptake(VFex) = -model.Vss(VFex);
if isfield(model,'Vuptake')
    Vuptake = model.Vuptake;
end

vh = [find(strcmpi(model.rxns,'exH'));...
      find(strcmpi(model.rxns,'exH2O'));...
      find(strcmpi(model.rxns,'exPI'))];

for irxn = 1:length(VFex)
%     sbid = logical(model.S(:,VFex(irxn))<0);
%     if ismember(VFex(irxn),vh)
%         ind = model.Vex(logical(model.S(sbid,model.Vex)));    
%         net_out = -sign(model.S(sbid,VFex(irxn)))*...
%                   (model.S(sbid,ind)*flux(ind));
%               -model.S(sbid,VFex(irxn))*Vuptake(VFex(irxn)));
        flux(VFex(irxn),:) = zeros(1,nc);%scale_flux(net_out);%-Vuptake(VFex(irxn)));
%         flux(VFex(irxn)) = scale_flux(flux(VFex(irxn)));
%     else
        
%     end
    
    
%     flux(VFex(irxn)) = scale_flux(kcatfwd(VFex(irxn))*mc(sbid)-...
%                        Vuptake(VFex(irxn)));    
%     vf = model.S(sbid,:);
%     vdiff = setdiff(find(vf),VFex(irxn));
%     if ~model.Vuptake(VFex(irxn))
%         flux(VFex(irxn)) = model.S(sbid,vdiff)*flux(vdiff);
%         %mc(sbid);
%     else
%         flux(VFex(irxn)) = -model.Vuptake(VFex(irxn));
%     end
end
flux = flux(VFex);
