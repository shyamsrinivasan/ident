function pvec = samplekcatRedox(model,pvec,mc)

R = 0.008314; %kJ/mol.K
T = 298; %K
RT = R*T;

% kcatfwd = pvec.kcat_fwd;
% kcatbkw = pvec.kcat_bkw;

vred = [find(strcmpi(model.rxns,'NADTRHD'))...
        find(strcmpi(model.rxns,'CYTBD'))...    %assumed kcat=100 s-1
        find(strcmpi(model.rxns,'THD2'))...     %assumed kcat=100 s-1
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
sigma = rand;
sb_rat = sigma/(1-sigma);
pvec.K(nadph,vred(1)) = 1e-5;%mc(nadph)/sb_rat;

% vflux(vred(1)) = kcatfwd(vred(1))*mc(nadph);

%CYTBD o2 --->
sigma = rand;
sb_rat = sigma/(1-sigma);
pvec.K(o2,vred(2)) = 1e-6;%mc(o2)/sb_rat;

% vflux(vred(2)) = kcatfwd(vred(2))*mc(o2);

%THD2 nadph ---> nadp
sigma = rand(2,1);
sb_rat = sigma(1)/(1-sigma(1));
pvec.K(nadph,vred(3)) = 1e-3;%mc(nadph)/sb_rat;
pr_rat = sigma(2)/(1-sigma(2));
pvec.K(nadp,vred(3)) = 1e-3;%mc(nadp)/pr_rat;

% vflux(vred(3)) = kcatfwd(vred(3))*mc(nadph) - kcatbkw(vred(3))*mc(nadp);

%NADH16 nad ---> nadh
sigma = rand(2,1);
sb_rat = sigma(1)/(1-sigma(1));
Ksb = mc(nad)/sb_rat;
Kscol(nad,1) = pvec.K(nad,vred(4));
if any(Kscol==1)
    pvec.K(Kscol==1,vred(4)) = 10;%Ksb(pvec.K(nad,vred(4))==1);
    pvec.Kind(Kscol==1,vred(4))=1;
end

pr_rat = sigma(2)/(1-sigma(2));
Kpr = mc(nadh)/pr_rat;
Kpcol(nadh,1) = pvec.K(nadh,vred(4));
if any(Kpcol==1)
    pvec.K(Kpcol==1,vred(4)) = 10;%Kpr(pvec.K(nadh,vred(4))==1);
    pvec.Kind(Kpcol==1,vred(4))=1;
end

if isnan(pvec.kcat_fwd(vred(4)))
    pvec.kcat_fwd(vred(4)) = pvec.kcat_bkw(vred(4))*model.Keq(vred(4))*...
                         pvec.K(nad,vred(4))/pvec.K(nadh,vred(4));
%     pvec.kcat_fwd(vred(4)) = pvec.kcat_bkw(vred(4))*exp(pvec.delGr(vred(4))/RT)*...
%                             (mc(nadh)/pvec.K(nadh,vred(4)))/...
%                             (mc(nad)/pvec.K(nad,vred(4)));
end
if isnan(pvec.kcat_bkw(vred(4)))
    pvec.kcat_bkw(vred(4)) = (pvec.kcat_fwd(vred(4))/model.Keq(vred(4)))*...
                              pvec.K(nadh,vred(4))/pvec.K(nad,vred(4));
end

% vflux(vred(4)) = kcatfwd(vred(4))*mc(nad) - kcatbkw(vred(4))*mc(nadh);

%ATPS4r atp ---> adp
sigma = rand(2,1);
sb_rat = sigma(1)/(1-sigma(1));
pvec.K(atp,vred(5)) = 1e-5;%mc(atp)/sb_rat;
pr_rat = sigma(2)/(1-sigma(2));
pvec.K(adp,vred(5)) = 1e-5;%mc(adp)/pr_rat;

if isnan(pvec.kcat_fwd(vred(5)))
    pvec.kcat_fwd(vred(5)) = pvec.kcat_bkw(vred(5))*model.Keq(vred(5))*...
                         pvec.K(atp,vred(5))/pvec.K(adp,vred(5));
end
if isnan(pvec.kcat_bkw(vred(5)))
    pvec.kcat_bkw(vred(5)) = (pvec.kcat_fwd(vred(5))/model.Keq(vred(5)))*...
                              pvec.K(adp,vred(5))/pvec.K(atp,vred(5));
end

if any(isnan(pvec.kcat_fwd(vred)))
    pvec.kcat_fwd(vred(isnan(pvec.kcat_fwd(vred)))) = 100;
end

if any(isnan(pvec.kcat_bkw(vred)))
   pvec.kcat_bkw(model.rev(vred(isnan(pvec.kcat_bkw(vred))))) = 100;
   pvec.kcat_bkw(~model.rev(vred(isnan(pvec.kcat_bkw(vred))))) = 0;   
end

pvec.Vmax(vred) = 1;
pvec.Vmax(model.Vss==0) = 0;

% vflux(vred(5)) = kcatfwd(vred(5))*mc(atp) - kcatbkw(vred(5))*mc(adp);