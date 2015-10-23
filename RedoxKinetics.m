function [flux,vflux,vred] = RedoxKinetics(model,pvec,mc,flux)
if nargin<4
    flux = zeros(model.nt_rxn,1);
end
vflux = zeros(length(model.Vss),1);

kcatfwd = pvec.kcat_fwd;
kcatbkw = pvec.kcat_bkw;
Vmax = pvec.Vmax;
K = pvec.K;

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

hc = strcmpi(model.mets,'h[c]');
he = strcmpi(model.mets,'h[e]');
q8h2 = strcmpi(model.mets,'q8h2[c]');
q8 = strcmpi(model.mets,'q8[c]');
pic = strcmpi(model.mets,'pi[c]');
h2o = strcmpi(model.mets,'h2o[c]');
    
%NADTRHD nadph --->
if mc(nad)>0 
    rat = mc(nadph)/K(nadph,vred(1));
    vflux(vred(1)) = kcatfwd(vred(1))*rat/(1+rat);
else
    vflux(vred(1)) = 0;
end
%CYTBD o2 --->
if mc(hc)>0 && mc(q8h2)>0
    rat = mc(o2)/K(o2,vred(2));
    vflux(vred(2)) = kcatfwd(vred(2))*rat/(1+rat);
else
    vflux(vred(2)) = 0;
end
%THD2 nadph ---> nadp
if mc(hc)>0 && mc(nad)>0
    ratf = mc(nadph)/K(nadph,vred(3));
    vfwd = kcatfwd(vred(3))*ratf;
%     vfwd = kcatfwd(vred(3))*ratf/(1+ratf);
else
    ratf = 0;
    vfwd = 0;
end
if mc(he)>0 && mc(nadh)>0
    ratb = mc(nadp)/K(nadp,vred(3));
    vbkw = kcatbkw(vred(3))*ratb;
%     vbkw = kcatbkw(vred(3))*ratb/(1+ratb);
else
    ratb = 0;
    vbkw = 0;
end
vflux(vred(3)) = (vfwd-vbkw)/(1+ratf+ratb);
% if mc(hc)>0 && mc(he)>0
%     vflux(vred(3)) = kcatfwd(vred(3))*mc(nadph)/K(nadph,vred(3)) -...
%                  kcatbkw(vred(3))*mc(nadp)/K(nadp,vred(3));
% elseif mc(hc)>0
%     vflux(vred(3)) = kcatfwd(vred(3))*mc(nadph)/K(nadph,vred(3));
% elseif mc(he)>0
%     vflux(vred(3)) = -kcatbkw(vred(3))*mc(nadp)/K(nadp,vred(3));
% end

%NADH16 nad ---> nadh
if mc(he)>0 && mc(q8h2)>0
    ratf = mc(nad)/K(nad,vred(4));
    vfwd1 = kcatfwd(vred(4))*ratf;%/(1+ratf);
else
    ratf = 0;
    vfwd1 = 0;
end
if mc(hc)>0 && mc(q8)>0
    ratb = mc(nadh)/K(nadh,vred(4));
    vbkw1 = kcatbkw(vred(4))*ratb;%/(1+ratb);
else
    ratb = 0;
    vbkw1 = 0;
end
vflux(vred(4)) = (vfwd1-vbkw1)/(1+ratf+ratb);
% vflux(vred(4)) = kcatfwd(vred(4))*mc(nad)/K(nad,vred(4)) -...
%                  kcatbkw(vred(4))*mc(nadh)/K(nadh,vred(4));

%ATPS4r atp ---> adp
%vfwd != 0 only under anaerobic conditions
if mc(hc)>0 && mc(h2o)>0
    ratf = mc(atp)/K(atp,vred(5));
    vfwd2 = kcatfwd(vred(5))*ratf;%/(1+ratf);
else
    ratf = 0;
    vfwd2 = 0;
end
vfwd2 = 0;
if mc(he)>0 && mc(pic)>0
    ratb = mc(adp)/K(adp,vred(5));
    vbkw2 = kcatbkw(vred(5))*ratb;%/(1+ratb);
else
    ratb = 0;
    vbkw2 = 0;
end
vflux(vred(5)) = (vfwd2-vbkw2)/(1+ratf+ratb);
% vflux(vred(5)) = kcatfwd(vred(5))*mc(atp)/K(atp,vred(5)) -...
%                  kcatbkw(vred(5))*mc(adp)/K(adp,vred(5));

flux(vred) = Vmax(vred).*scale_flux(vflux(vred));
% flux = flux(vred);
vflux = vflux(vred);