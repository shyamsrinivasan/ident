function [flux,vflux] = CKinetics(model,pvec,mc,Vind)
S = model.S;
nrxn = model.nt_rxn;
K = pvec.K;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

vflux = zeros(nrxn,1); 
flux = zeros(nrxn,1);

%eliminate consideration for excess cofators
%pi[c],pi[e],h[c],h[e],h2o[c]
he = find(strcmpi(model.mets,'h[e]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
vmet = [he...
        find(strcmpi(model.mets,'h[c]'))...
        find(strcmpi(model.mets,'pi[e]'))...
        find(strcmpi(model.mets,'pi[c]'))...
        h2o...       
        find(strcmpi(model.mets,'co2[c]'))];


for irxn = 1:length(Vind)
    %compensated species indices
%     sbcmp = zeros(length(model.mets),1);
%     prcmp = zeros(length(model.mets),1);
    nr_flux = zeros(1,1);
    
    nmet = size(S,1);
    
    sbid = S(:,Vind(irxn))<0;
    prid = S(:,Vind(irxn))>0;    
    
    if ~any(strcmpi(model.rxns{Vind(irxn)},'PPC'))
        if any(sbid)
            mc_alls = prod(logical(mc(sbid)));
            if ~any(model.CMPS(sbid,Vind(irxn)))            
                sbid(vmet) = 0;
                cmp_s = [];
            else
                sbid = find(sbid);
                cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
                sbid = setdiff(sbid,cmp_s);
                sbid = setdiff(sbid,[he h2o]);
                sbid = logical(sparse(sbid,1,1,nmet,1));
    %             sbid(sbid==cmp_s) = [];                     
            end
        end

        if any(prid)
            mc_allp = prod(logical(mc(prid)));
            if ~any(model.CMPS(prid,Vind(irxn)))
                prid(vmet) = 0;
                cmp_p = [];
            else
                prid = find(prid);
                cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
                prid = setdiff(prid,cmp_p);  
                prid = setdiff(prid,[he h2o]);
                prid = logical(sparse(prid,1,1,nmet,1));
            end
        end        
    else
        mc_alls = prod(logical(mc(sbid)));
        mc_allp = prod(logical(mc(prid)));
        sbid(h2o) = 0;
        prid(h2o) = 0;        
        if ~any(model.CMPS(sbid,Vind(irxn)))  
            sbid(vmet) = 0;
            cmp_s = [];
        else
            sbid = find(sbid);
            cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
            sbid = setdiff(sbid,cmp_s);
            sbid = setdiff(sbid,[he h2o]);    
            sbid = logical(sparse(sbid,1,1,nmet,1));
        end
        if ~any(model.CMPS(prid,Vind(irxn)))
            prid(vmet) = 0;
            cmp_p = [];
        else
            prid = find(prid);
            cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
            prid = setdiff(prid,cmp_p);  
            prid = setdiff(prid,[he h2o]);
            prid = logical(sparse(prid,1,1,nmet,1));
        end             
    end
    if ~isempty(cmp_s)
        cmp_s = prod(mc(cmp_s).*(-model.S(cmp_s,Vind(irxn))));
%         if cmp_s > 0
%             cmp_s = 1;
%         else
%             cmp_s = 0;
%         end
    else
        cmp_s = 1;
    end
    if ~isempty(cmp_p)
        cmp_p = prod(mc(cmp_p).*(model.S(cmp_p,Vind(irxn))));
%         if cmp_p > 0
%             cmp_p = 1;
%         else
%             cmp_p = 0;
%         end
    else
        cmp_p = 1;
    end
%         if any(prid(vmet))
%             prcmp(vmet(logical(prid(vmet)))) =...
%             prid(vmet(logical(prid(vmet))));
%             prid(vmet) = 0;
%         end
%     end
    
%     cmp_s = model.CMPS(sbid,Vind(irxn));
%     cmp_s = prod(mc(model.CMPS(sbid,Vind(irxn))).^model.S(
%     cmp_s = prod(mc(logical(sbcmp)).^sbcmp(logical(sbcmp)));
%     cmp_p = prod(mc(logical(prcmp)).^prcmp(logical(prcmp)));
    
    if model.rev(Vind(irxn))
        Sb = -S(sbid,Vind(irxn));
        Sp = S(prid,Vind(irxn));   
        if all(mc(sbid)>0) && all(mc(prid)>0)
            nr_flux = mc_alls*kfwd(Vind(irxn))*cmp_s*prod(mc(sbid)./K(sbid,Vind(irxn))) -...
                      mc_allp*kbkw(Vind(irxn))*cmp_p*prod(mc(prid)./K(prid,Vind(irxn)));
        elseif all(mc(sbid)>0)
            nr_flux = mc_alls*kfwd(Vind(irxn))*cmp_s*prod(mc(sbid)./K(sbid,Vind(irxn)));
        elseif all(mc(prid)>0)
            nr_flux = -mc_allp*kbkw(Vind(irxn))*cmp_p*prod(mc(prid)./K(prid,Vind(irxn)));
        end            
    elseif ~model.rev(Vind(irxn))
        Sb = -S(sbid,Vind(irxn));    
        Sp = S(prid,Vind(irxn));
        if all(mc(sbid)>0)
            nr_flux =...
            mc_alls*kfwd(Vind(irxn))*cmp_s*prod(mc(sbid)./K(sbid,Vind(irxn)));
        end
    end
    
    if any(sbid) && any(prid)
        
        %Numerator - 1.6
%         if model.rev(Vind(irxn))%reversible
%             if all(mc(sbid)>=0)&&all(mc(prid)>=0)
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 prid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));
%                 Sp = -S(prid,Vind(irxn));
%                 nr_flux = kfwd(Vind(irxn))*prod(mc(sbid)./K(sbid,Vind(irxn))) -...
%                           kbkw(Vind(irxn))*prod(mc(prid)./K(prid,Vind(irxn)));
%             elseif all(mc(sbid)>=0)
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 prid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));
%                 Sp = -S(prid,Vind(irxn));
%                 %Numerator - 1.6
%                 nr_flux = kfwd(Vind(irxn))*prod(mc(sbid)./K(sbid,Vind(irxn)));
%             elseif all(mc(prid)>=0)
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 prid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));
%                 Sp = -S(prid,Vind(irxn));
%                 nr_flux = -kbkw(Vind(irxn))*prod(mc(prid)./K(prid,Vind(irxn)));
%             else
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 prid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));
%                 Sp = -S(prid,Vind(irxn));
%                 nr_flux = 0;
%             end
%         elseif ~model.rev(Vind(irxn))%irreversible
%             if all(mc(sbid)>=0)
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));    
%                 Sp = -S(prid,Vind(irxn));
%                 nr_flux = kfwd(Vind(irxn))*prod(mc(sbid)./K(sbid,Vind(irxn)));
%             else
%                 sbid([pic pie hc he h2o co2]) = 0;
%                 Sb = S(sbid,Vind(irxn));
%                 Sp = -S(prid,Vind(irxn));
%                 nr_flux = 0;
%             end
%         end
        
        %Denominator - 1.6
        dr_sb = 1+mc(sbid)./K(sbid,Vind(irxn));
        for j = 1:length(find(sbid))
            for si = 2:Sb(j)
                dr_sb(j) = dr_sb(j) + dr_sb(j)^si;                
            end
        end
        %dr_pr
        dr_pr = 1+mc(prid)./K(prid,Vind(irxn));
        for j = 1:length(find(prid))
            for si = 2:Sp(j)
                dr_pr(j) = dr_pr(j)+dr_pr(j)^si;
            end
        end
        dr_flux = prod(dr_sb)+prod(dr_pr)-1;
        vflux(Vind(irxn)) = scale_flux(nr_flux/dr_flux);
    else
        vflux(Vind(irxn)) = 0;
    end
    flux(Vind(irxn)) = Vmax(Vind(irxn))*vflux(Vind(irxn));
end
flux = flux(Vind);
vflux = vflux(Vind);