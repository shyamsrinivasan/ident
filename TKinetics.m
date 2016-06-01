function [flux,vflux] = TKinetics(model,pvec,mc,Vex)
flux = zeros(model.nt_rxn,1);
vflux = zeros(model.nt_rxn,1);
Vmax = pvec.Vmax;
K = pvec.K;
S = model.S;

he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
% h2o = find(strcmpi(model.mets,'h2o[c]'));
% pie = find(strcmpi(model.mets,'pi[e]'));
% pic = find(strcmpi(model.mets,'pi[c]'));
% co2 = find(strcmpi(model.mets,'co2[c]'));
    
% q8 = find(strcmpi(model.mets,'q8[c]'));
% q8h2 = find(strcmpi(model.mets,'q8h2[c]'));

% vmet = [he hc pie pic h2o co2];   

% vpts = find(strcmpi(model.rxns,'GLCpts'));
% vthd2 = find(strcmpi(model.rxns,'THD2'));%nadph --->
% vatps = find(strcmpi(model.rxns,'ATPS4r'));%atp ---> adp
% vnad16 = find(strcmpi(model.rxns,'NADH16'));%nad ---> nadh
% vcyt = find(strcmpi(model.rxns,'CYTBD'));%o2 --->
% vspl = [vpts vthd2 vatps vnad16 vcyt];

% nadph = strcmpi(model.mets,'nadph[c]');
% nadp = strcmpi(model.mets,'nadp[c]');
% nad = strcmpi(model.mets,'nad[c]');
% nadh = strcmpi(model.mets,'nadh[c]');
% atp = strcmpi(model.mets,'atp[c]');
% adp = strcmpi(model.mets,'adp[c]');
% o2 = strcmpi(model.mets,'o2[c]');

for irxn = 1:length(Vex)    
        
    %kcat
    kfwd = pvec.kcat_fwd(Vex(irxn));
    kbkw = pvec.kcat_bkw(Vex(irxn));

    nmet = size(S,1);
    %kinetics - substrate and product
    sbid = logical(model.S(:,Vex(irxn))<0);
    prid = logical(model.S(:,Vex(irxn))>0); 
    
    % remove protons
    sbid([he hc]) = 0;
    prid([he hc]) = 0;
    
    if any(strcmpi(model.rxns{Vex(irxn)},'O2t'))
        nr_flux = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
                  (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
                  prod(mc(prid)./K(prid,Vex(irxn))));
    elseif any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
        Kapie = 0.89; % mM
        nr_flux = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
                  (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
                  prod(mc(prid)./K(prid,Vex(irxn))))*...
                  1/(1+Kapie/mc(sbid));
    else
        if model.rev(Vex(irxn))
            nr_flux = kfwd*(prod(mc(sbid)./K(sbid,Vex(irxn)))-...
                      prod(mc(prid)./K(prid,Vex(irxn))))/...
                      (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
                      prod(mc(prid)./K(prid,Vex(irxn))));
        elseif ~model.rev(Vex(irxn))
            nr_flux = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
                      (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
                      prod(mc(prid)./K(prid,Vex(irxn))));
        end
    end
    if any(sbid) && any(prid)
        vflux(Vex(irxn)) = scale_flux(nr_flux);
    else
        vflux(Vex(irxn)) = 0;
    end
 
%     if model.rev(Vex(irxn))         
%         if all(mc(sbid)>0) && all(mc(prid)>0)
%             nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn))) -...
%                       mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
%         elseif all(mc(sbid)>0)
%             nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));
%         elseif all(mc(prid)>0)
%             nr_flux = -mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
%         end
%         if any(sbid) && any(prid)
%             % Denominator - 1.6
%             dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));            
%             % dr_pr
%             dr_pr = 1+mc(prid)./K(prid,Vex(irxn));            
%         else
%             dr_sb = 0;
%             dr_pr = 0;
%         end        
%     elseif ~model.rev(Vex(irxn))        
%         nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));  
%         if any(sbid) && any(prid)
%             dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));            
%         else
%             dr_sb = 0;
%         end
%     end
%     
%     if any(sbid) && any(prid)
%         %Denominator - 1.6
%         dr_flux = prod(dr_sb)+prod(dr_pr)-1;
%         vflux(Vex(irxn)) = scale_flux(nr_flux/dr_flux);
%     else
%         vflux(Vex(irxn)) = 0;
%     end
    
%     if any(strcmpi(model.rxns{Vex(irxn)},'o2t'))
% %         Do2 = 2.1e-9; %m2/s
% %         Acell = 4.42e-12; %m2
% %         vflux(Vex(irxn)) = Do2/Acell*(mc(sbid)-mc(prid));
%         vflux(Vex(irxn)) = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
%                            (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
%                            prod(mc(prid)./K(prid,Vex(irxn))));
% %         Vmax(Vex(irxn)) = 1;
%     end
%     flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
    vflux(Vex(irxn)) = scale_flux(vflux(Vex(irxn)));
    flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
end
flux = flux(Vex);
vflux = vflux(Vex);
