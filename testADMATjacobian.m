%function for example using ADMAT
function [y,J] = testADMATjacobian(data,mc)
addpath(genpath('C:\Users\shyam\Documents\MATLAB\zz_ADMAT-2.0'),'-end');
mc_obj = deriv(mc,eye(data.nt_metab));
fx = rhsfunc(mc_obj,data);
y = getval(fx);
J = getydot(fx);
rmpath('C:\Users\shyam\Documents\MATLAB\zz_ADMAT-2.0');
return

function rhs = rhsfunc(x,data)
    pvec = data.pvec;
    nt_m = data.nt_metab;
    nin_m = data.nint_metab;
    rhs = zeros(nt_m,1);%[Metabolites;Biomass]
    rhs = cons(rhs,x);
    flux = allflux(data,pvec,x);

    %calculate f(x,p) in dXdt = f(x,p)
    rhs(1:nin_m) = data.S(1:nin_m,:)*flux;
    rhs(nin_m+1:nt_m) = data.S(nin_m+1:nt_m,:)*zeros(data.nt_rxn,1);
return

function flux = allflux(data,pvec,x)
    flux = zeros(data.nt_rxn,1);
    flux = cons(flux,x);
    Vind = data.Vind;
    Vex = setdiff(data.Vex,Vind);
    VFex = data.VFex;
    flux(Vind) = CKin(data,pvec,x,Vind);
    flux(Vex) = TKin(data,pvec,x,Vex);
    flux(VFex) = EKin(data,pvec,x,VFex,flux);
    flux(strcmpi(data.rxns,'atpm')) =...
    8.39*logical(getval(x(strcmpi(data.mets,'atp[c]'))));
    flux(data.bmrxn) = data.Vss(data.bmrxn);
return

function flux = CKin(model,pvec,mc,Vind)
nrxn = model.nt_rxn;
vflux = zeros(nrxn,1); 
vflux =cons(vflux,mc);
flux = zeros(nrxn,1);
flux = cons(flux,mc);

S = model.S;
K = pvec.K;
kfwd = pvec.kcat_fwd;
kbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

for irxn = 1:length(Vind)
    nr_flux = zeros(1,1);
    nr_flux = cons(nr_flux,mc);
    nmet = size(S,1);
    
    sbid = S(:,Vind(irxn))<0;
    prid = S(:,Vind(irxn))>0;
    
    if any(sbid)
        mc_alls = prod(logical(getval(mc(sbid))));
        if ~any(model.CMPS(sbid,Vind(irxn)))            
            cmp_s = [];
        else
            sbid = find(sbid);
            cmp_s = sbid(logical(model.CMPS(sbid,Vind(irxn))));
            sbid = setdiff(sbid,cmp_s);            
            sbid = logical(sparse(sbid,1,1,nmet,1));                
        end
    end

    if any(prid)
        mc_allp = prod(logical(getval(mc(prid))));
        if ~any(model.CMPS(prid,Vind(irxn)))            
            cmp_p = [];
        else
            prid = find(prid);
            cmp_p = prid(logical(model.CMPS(prid,Vind(irxn))));
            prid = setdiff(prid,cmp_p);  
            prid = logical(sparse(prid,1,1,nmet,1));
        end
    end   
    if ~isempty(cmp_s)
        cmp_s = prod(mc(cmp_s).*(-model.S(cmp_s,Vind(irxn))));
    else
        cmp_s = 1;
    end
    if ~isempty(cmp_p)
        cmp_p = prod(mc(cmp_p).*(model.S(cmp_p,Vind(irxn))));
    else
        cmp_p = 1;
    end
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
    if any(sbid)&&any(prid)
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
        vflux(Vind(irxn)) = scaleF(nr_flux/dr_flux);
    else
        vflux(Vind(irxn)) = 0;
    end
    %final flux(Vind(irxn))
    flux(Vind(irxn)) = Vmax(Vind(irxn))*vflux(Vind(irxn));        
end
flux = flux(Vind);
return

function flux = TKin(model,pvec,mc,Vex)
nrxn = model.nt_rxn;
vflux = zeros(nrxn,1); 
vflux =cons(vflux,mc);
flux = zeros(nrxn,1);
flux = cons(flux,mc);

Vmax = pvec.Vmax;
K = pvec.K;
S = model.S;

for irxn = 1:length(Vex)
    kfwd = pvec.kcat_fwd(Vex(irxn));
    kbkw = pvec.kcat_bkw(Vex(irxn));
    
    nr_flux = zeros(1,1);
    nr_flux = cons(nr_flux,mc);
    
    nmet = size(S,1);    
    sbid = logical(model.S(:,Vex(irxn))<0);
    prid = logical(model.S(:,Vex(irxn))>0);
    
    if any(sbid)
        mc_alls = prod(logical(getval(mc(sbid))));              
        if ~any(model.CMPS(sbid,Vex(irxn)))    
            cmp_s = [];
        else
            sbid = find(sbid);
            cmp_s = sbid(logical(model.CMPS(sbid,Vex(irxn))));
            sbid = setdiff(sbid,cmp_s);
            sbid = logical(sparse(sbid,1,1,nmet,1));
        end
    end
    if any(prid)
        mc_allp = prod(logical(getval(mc(prid))));        
        if ~any(model.CMPS(prid,Vex(irxn)))
            cmp_p = [];
        else
            prid = find(prid);
            cmp_p = prid(logical(model.CMPS(prid,Vex(irxn))));
            prid = setdiff(prid,cmp_p);
            prid = logical(sparse(prid,1,1,nmet,1));
        end
    end 
    if ~isempty(cmp_s)
        cmp_s = prod(mc(cmp_s).*(-model.S(cmp_s,Vex(irxn))));
    else
        cmp_s = 1;
    end
    if ~isempty(cmp_p)
        cmp_p = prod(mc(cmp_p).*(model.S(cmp_p,Vex(irxn))));
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
    elseif ~model.rev(Vex(irxn))
        Sb = -S(sbid,Vex(irxn));    
        Sp = S(prid,Vex(irxn));
        if all(mc(sbid)>0)
            nr_flux = mc_alls*kfwd*cmp_s*prod(mc(sbid)./K(sbid,Vex(irxn)));
        end
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
        dr_flux = prod(dr_sb)+prod(dr_pr)-1;
        vflux(Vex(irxn)) = scaleF(nr_flux/dr_flux);
    else
        vflux(Vex(irxn)) = 0;
    end
    vflux(Vex(irxn)) = scaleF(vflux(Vex(irxn)));
    flux(Vex(irxn)) = Vmax(Vex(irxn))*vflux(Vex(irxn));
end
flux = flux(Vex);
return

function flux = EKin(model,pvec,mc,VFex,flux)
%this function is included just to be consistent with what I am doing on my
%end
for irxn = 1:length(VFex)
    flux(VFex(irxn)) = 0;
end
flux = flux(VFex);
return

function flux = scaleF(flux)
    for irxn = 1:length(flux)
        vf = flux(irxn);
        if ~isnan(vf) && getval(vf)
            if abs(vf)<=1e-10
                flux(irxn) = 0;
            end    
        end
    end
return
