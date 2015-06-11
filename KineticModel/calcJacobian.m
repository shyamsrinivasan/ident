function [Jsat,Jtherm,Jreg,...
          relJsat,relJtherm,relJreg] = calcJacobian(subs,pruds,act,inhib,...
                                                      Ksubs,Kprod,Kact,Kinb,...
                                                      s_subs,s_prod,s_act,s_inhib,...
                                                      numrsat,drsatsubs,drsatprod,...
                                                      vact,vinhib,...
                                                      gamma)
% nrxn = model.nrxn;
% nmetab = model.nmetab;

%Inputs: Ksubs, Kprod, KIact, KIinb, s_subs, s_prod, s_act, s_inhib,
%        subs, pruds
%Outputs: jacobians & relative jacobinas
%i.e.
%Jacobians: dv/dS = Jsat+Jtherm+Jreg
%           dv/dS = v0/S0*(relJsat+relJtherm+relJreg) - Not complete yet
   
%Saturation Component
subsvec1 = s_subs.*(subs.^(s_subs - 1))./(subs.^s_subs);
subsvec2 = (s_subs./Ksubs).*(((1 + subs./Ksubs).^(s_subs-1))./((1 + subs./Ksubs).^s_subs));
subsscalar = drsatsubs*drsatprod/(drsatsubs*drsatprod - 1);
Jsat_subs = (numrsat/(drsatsubs*drsatprod - 1))*(subsvec1 - subsscalar*subsvec2);
relJsat_subs = subs.*(subsvec1 - subsscalar*subsvec2);

prodscalar = drsatsubs*drsatprod;
prodvector = (s_prod./Kprod).*(((1 + pruds./Kprod).^(s_prod-1))./((1 + pruds./Kprod).^s_prod));
Jsat_prod = -(numrsat/(drsatsubs*drsatprod-1)^2)*prodscalar*prodvector;
relJsat_prod = -(pruds/(drsatsubs*drsatprod - 1)).*(prodscalar*prodvector);

Jsat = [Jsat_subs' Jsat_prod'];         %jacobians - dv/dS
relJsat = [relJsat_subs' relJsat_prod'];%relative jacobians - dlnv/dlnS

%Thermodynamic Component
subsvec = (s_subs.*(subs.^(s_subs - 1)))./(subs.^s_subs);
Jtherm_subs = gamma*subsvec;
relJtherm_subs = (gamma/(1-gamma))*(subs.*subsvec);

prodvec = (s_prod.*(pruds.^(s_prod - 1)))./(pruds.^s_prod);
Jtherm_prod = - gamma*prodvec;
relJtherm_prod = (gamma/(1-gamma))*(pruds.*prodvec);

Jtherm = [Jtherm_subs' Jtherm_prod'];
relJtherm = [relJtherm_subs' relJtherm_prod'];

%Enzyme Regulation
act_ratio = (act.^s_act)./(Kact.^s_act);
Jreg_act = ((s_act./act).*(1-act_ratio./(1+act_ratio)));
%relJreg_act = act.*Jreg_act;

Jreg_inhib = -(s_inhib.*(inhib.^s_inhib-1)./(Kinb.^s_inhib));
%relJreg_inhib = inhib.*Jreg_inhib;

Jreg = vact*vinhib*[Jreg_act' Jreg_inhib'];
relJreg = [(act.*Jreg_act)' (inhib.*Jreg_inhib)'];

%J = Jsat + Jtherm + Jreg;   


%end
