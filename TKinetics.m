function [flux,vflux] = TKinetics(model,pvec,mc,Vex)
flux = zeros(model.nt_rxn,1);
vflux = zeros(model.nt_rxn,1);
Vmax = pvec.Vmax;
K = pvec.K;
S = model.S;

he = find(strcmpi(model.mets,'h[e]'));
hc = find(strcmpi(model.mets,'h[c]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
pie = find(strcmpi(model.mets,'pi[e]'));
pic = find(strcmpi(model.mets,'pi[c]'));
co2 = find(strcmpi(model.mets,'co2[c]'));
    
q8 = find(strcmpi(model.mets,'q8[c]'));
q8h2 = find(strcmpi(model.mets,'q8h2[c]'));

vmet = [he hc pie pic h2o co2];   

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
%     if ~isnan(pvec.kcat_fwd(Vex(irxn)))
%         kfwd = pvec.kcat_fwd(Vex(irxn));
%     else
%         kfwd = 1000;
%     end
%     if ~isnan(pvec.kcat_bkw(Vex(irxn)))
%         kbkw = pvec.kcat_bkw(Vex(irxn));
%     else
%         kbkw = 1000;
%     end 
    nmet = size(S,1);
    %kinetics - substrate and product
    sbid = logical(model.S(:,Vex(irxn))<0);
    prid = logical(model.S(:,Vex(irxn))>0); 
    
    if any(sbid)
        mc_alls = prod(logical(mc(sbid)));
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot'))
%             sbid(h2o) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             sbid([pie pic]) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'CO2t'))
%             sbid(co2) = 0;
%         end
%         
        if ~any(model.CMPS(sbid,Vex(irxn))) 
            sbid([he hc]) = 0;     
            cmp_s = [];
        else
            sbid = find(sbid);
            cmp_s = sbid(logical(model.CMPS(sbid,Vex(irxn))));
            sbid = setdiff(sbid,cmp_s);
            sbid = logical(sparse(sbid,1,1,nmet,1));
%             sbid = setdiff(sbid,[he h2o]);
        end
    end    

    if any(prid)
        mc_allp = prod(logical(mc(prid)));
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot'))
%             prid(h2o) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             prid([pie pic]) = 0;
%         end
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'CO2t'))
%             prid(co2) = 0;
%         end
        if ~any(model.CMPS(prid,Vex(irxn)))
            prid([he hc]) = 0;
            cmp_p = [];
        else
            prid = find(prid);
            cmp_p = prid(logical(model.CMPS(prid,Vex(irxn))));
            prid = setdiff(prid,cmp_p);
            prid = logical(sparse(prid,1,1,nmet,1));
%             prid = setdiff(prid,[he h2o]);
        end
    end 

    if ~isempty(cmp_s)
        cmp_s = logical(prod(mc(cmp_s).*(-model.S(cmp_s,Vex(irxn)))));
    else
        cmp_s = 1;
    end
    if ~isempty(cmp_p)
        cmp_p = logical(prod(mc(cmp_p).*(model.S(cmp_p,Vex(irxn)))));
    else
        cmp_p = 1;
    end
    
    if model.rev(Vex(irxn))
        Sb = -S(sbid,Vex(irxn));
        Sp = S(prid,Vex(irxn));  
        if all(mc(sbid)>0) && all(mc(prid)>0)
            nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn))) -...
                      mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
        elseif all(mc(sbid)>0)
            nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));
        elseif all(mc(prid)>0)
            nr_flux = -mc_allp*kbkw*cmp_p*prod(mc(prid)./K(prid,Vex(irxn)));
        end
        if any(sbid) && any(prid)
            %Denominator - 1.6
            dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));
            for j = 1:length(find(sbid))
                for si = 2:Sb(j)
                    dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                end
            end
            %dr_pr
            dr_pr = 1+mc(prid)./K(prid,Vex(irxn));
            for j = 1:length(find(prid))
                for si = 2:Sp(j)
                    dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
                end
            end
        else
            dr_sb = 0;
            dr_pr = 0;
        end        
    elseif ~model.rev(Vex(irxn))
        Sb = -S(sbid,Vex(irxn));    
        Sp = S(prid,Vex(irxn));
        nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));  
        if any(sbid) && any(prid)
            dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));
            for j = 1:length(find(sbid))
                for si = 2:Sb(j)
                    dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
                end
            end
        else
            dr_sb = 0;
        end
    end
    
    if any(sbid) && any(prid)
        %Denominator - 1.6
%         dr_sb = 1+mc(sbid)./K(sbid,Vex(irxn));
%         for j = 1:length(find(sbid))
%             for si = 2:Sb(j)
%                 dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
%             end
%         end
        %dr_pr
%         dr_pr = 1+mc(prid)./K(prid,Vex(irxn));
%         for j = 1:length(find(prid))
%             for si = 2:Sp(j)
%                 dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
%             end
%         end
        dr_flux = prod(dr_sb)+prod(dr_pr)-1;
        vflux(Vex(irxn)) = scale_flux(nr_flux/dr_flux);
    else
        vflux(Vex(irxn)) = 0;
    end
    
    if any(strcmpi(model.rxns{Vex(irxn)},'o2t'))
%         Do2 = 2.1e-9; %m2/s
%         Acell = 4.42e-12; %m2
%         vflux(Vex(irxn)) = Do2/Acell*(mc(sbid)-mc(prid));
        vflux(Vex(irxn)) = kfwd*prod(mc(sbid)./K(sbid,Vex(irxn)))/...
                           (1+prod(mc(sbid)./K(sbid,Vex(irxn)))+...
                           prod(mc(prid)./K(prid,Vex(irxn))));
%         Vmax(Vex(irxn)) = 1;
    end
%     flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));

        
%     if model.rev(Vex(irxn))
%         vfwd = mc_alls*cmp_s*kfwd*sr_rt;
%         vbkw = mc_allp*cmp_p*kbkw*pr_rt;
%         vflux(Vex(irxn)) = (vfwd-vbkw)/(1+sr_rt+pr_rt);
%     elseif ~model.rev(Vex(irxn))
%         vfwd = mc_alls*cmp_s*kfwd*sr_rt;
%         vflux(Vex(irxn)) = vfwd/(1+sr_rt);
%     end
    
%     if all(mc(sbid)>0)
%         %all forward 
%         %eliminate water from all reactions
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot')) &&...
%            ~any(strcmpi(model.rxns{Vex(irxn)},'ATPS4r')) &&...
%            ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))            
%             sbid([h2o pic pie hc he q8 q8h2]) = 0;
%         elseif ~any(strcmpi(model.rxns{Vex(irxn)},'ATPS4r')) &&...
%                ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%            sbid([pic pie hc he q8 q8h2]) = 0;
%         elseif ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             sbid([h2o pie hc he q8 q8h2]) = 0;
%         else
%             sbid([h2o hc he q8 q8h2]) = 0;
%         end            
%         sr_rt=min(mc(sbid)./K(sbid,Vex(irxn)));            
%     else
%         %all foward = 0
%         sr_rt = 0;
%     end
%     vfwd = kfwd*min(sr_rt./(sr_rt+1));
    
%     if all(mc(prid)>0)
%         %all reverse
%         %eliminate water from all reactions
%         if ~any(strcmpi(model.rxns{Vex(irxn)},'h2ot')) &&...
%            ~any(strcmpi(model.rxns{Vex(irxn)},'ATPS4r')) &&...
%            ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))            
%             prid([h2o pic pie hc he q8 q8h2]) = 0;
%         elseif ~any(strcmpi(model.rxns{Vex(irxn)},'ATPS4r')) &&...
%                ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             prid([pic pie hc he q8 q8h2]) = 0;
%         elseif ~any(strcmpi(model.rxns{Vex(irxn)},'PIt2r'))
%             prid([h2o pie hc he q8 q8h2]) = 0;
%         else
%             prid([h2o hc he q8 q8h2]) = 0;
%         end     
%         pr_rt = min(mc(prid)./K(prid,Vex(irxn)));
%     else
%         %all reverse = 0
%         pr_rt = 0;
%     end
%     vbkw = kbkw*min(pr_rt./(pr_rt+1));
       
%     if ~ismember(Vex(irxn),vspl)
%         if all(mc(sbid))
%             sr_rt=mc(sbid)./(0.9*mc(sbid));
%         else
%             sr_rt=mc(sbid)./9e-9;
%         end
% 
%         if all(mc(prid))
%             pr_rt = mc(prid)./(0.9*mc(prid));
%         else
%             pr_rt=mc(prid)./9e-9;
%         end
        
%         if model.rev(Vex(irxn))                
% %             vflux(Vex(irxn)) = kfwd*min(sr_rt./(sr_rt+1)) -...
% %                                kbkw*min(pr_rt./(pr_rt+1));
%             vflux(Vex(irxn)) = (kfwd*sr_rt-kbkw*pr_rt)/(1+sr_rt+pr_rt);
% %             vflux(Vex(irxn)) = vfwd-vbkw;
%                            
%         elseif ~model.rev(Vex(irxn))
% %             vflux(Vex(irxn)) =...
% %             kfwd*min(sr_rt./(sr_rt+1));
%             vflux(Vex(irxn)) = kfwd*sr_rt/(1+sr_rt);
% %             vflux(Vex(irxn)) = vfwd;
%         end  
%         fprintf('%d. Vmax = %3.6g \tvflux = %3.6g\n',...
%                 Vex(irxn),Vmax(Vex(irxn)),vflux(Vex(irxn)));
%     end
    
    %other fluxes
%     if Vex(irxn)==vpts
%         [~,vflux(Vex(irxn))] = CKinetics(model,pvec,mc,Vex(irxn));
%     end
%     if Vex(irxn)==vthd2
%         vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
%         %kfwd*mc(nadph)-kbkw*mc(nadp);
%     end
%     if Vex(irxn)==vatps
%         vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
%         %kfwd*mc(atp)-kbkw*mc(adp);
%     end
%     if Vex(irxn)==vnad16
%         vflux(Vex(irxn)) = CKinetics(model,pvec,mc,Vex(irxn));
%         %kfwd*mc(nad)-kbkw*mc(nadh);
%     end
%     if Vex(irxn)==vcyt
%         vflux(Vex(irxn)) = kfwd*(mc(o2))^0.5;
%     end   
    vflux(Vex(irxn)) = scale_flux(vflux(Vex(irxn)));
    flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
end
flux = flux(Vex);
vflux = vflux(Vex);
