function [flux,vflux] = RedoxKinetics(model,pvec,mc)

kcatfwd = pvec.kcat_fwd;
kcatbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;

%other reactions - redox balance
vred = [find(strcmpi(model.rxns,'NADTRHD'))...
        find(strcmpi(model.rxns,'CYTBD'))...
        find(strcmpi(model.rxns,'THD2'))...
        find(strcmpi(model.rxns,'NADH16'))...
        find(strcmpi(model.rxns,'ATPS4r'))];
    
nadph = strcmpi(model.mets,'nadph[c]');    
o2 = strcmpi(model.mets,'o2[c]');
nadp = strcmpi(model.mets,'nadp[c]');
nad = strcmpi(model.mets,'nad[c]');
nadh = strcmpi(model.mets,'nadh[c]');
atp = strcmpi(model.mets,'atp[c]');
adp = strcmpi(model.mets,'adp[c]');
    
%NADTRHD nadph --->
vflux(vred(1)) = kcatfwd(vred(1))*mc(nadph);
%CYTBD o2 --->
vflux(vred(2)) = kcatfwd(vred(2))*mc(o2);
%THD2 nadph ---> nadp
vflux(vred(3)) = kcatfwd(vred(3))*mc(nadph) - kcatbkw(vred(3))*mc(nadp);
%NADH16 nad ---> nadh
vflux(vred(4)) = kcatfwd(vred(4))*mc(nad) - kcatbkw(vred(4))*mc(nadh);
%ATPS4r atp ---> adp
vflux(vred(5)) = kcatfwd(vred(5))*mc(atp) - kcatbkw(vred(5))*mc(adp);

flux(vred) = Vmax(vred).*vflux(vred);