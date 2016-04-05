function flux = ETCflux(model,mc,flux)
if nargin < 3
    flux = zeros(length(model.rxns),1);
end
mc = mc.*1e-3;
ETCrxn = {'ATPS4r','NADH16','CYTBD','SUCDi','FRD7'};
hc = find(strcmpi(model.mets,'h[c]'));
he = find(strcmpi(model.mets,'h[e]'));
h2o = find(strcmpi(model.mets,'h2o[c]'));
q8 = strcmpi(model.mets,'q8[c]');
q8h2 = find(strcmpi(model.mets,'q8h2[c]'));
%deltapH
delpH = -log10(mc(hc))+log10(mc(he));%-deltapH
%pmf
psi_con = 0.7;%0-1
R = 0.008314;%kJ/mol.K
T = 298;%K
F = 0.0965;%kJ/mol.mV
Z = R*T/F;
pmf = (2.303*R*T/F)*delpH*1/(1-psi_con);
delG_atphydro = -30.5;%kJ/mol.K
%redox potentials in mV
Eq8 = 70;
Enad_nadh = -320;
EUbiOxi = 250; 

for i = 1:length(ETCrxn)
    tfr = strcmpi(model.rxns,ETCrxn{i});
    if any(tfr)
        sbid = model.S(:,tfr)<0;
        prid = model.S(:,tfr)>0;
        sbid([hc he h2o]) = 0;
        prid([hc he h2o]) = 0;
        switch(model.rxns{tfr})
            case 'ATPS4r'
                % adp + pi + 3he + hc <---> atp + h2o + 3hc
                gamma = 4*pmf*F+delG_atphydro-R*T*log(prod(mc(sbid))/prod(mc(prid)));
                Katps = 20; %mmol/gDCW.h
                flux(tfr) = -Katps/3600*(exp(gamma/(R*T))-1)/(exp(gamma/(R*T))+1);
            case 'NADH16'
                % nadh + q8 + 5hc <---> nad + q8h2 + 4he
                Knadh16 = 3.35; %mmole/gDCW.h.mV
                flux(tfr) = Knadh16/3600*(Eq8-Enad_nadh-Z/2*log(prod(mc(sbid))/prod(mc(prid)))-2*pmf);
            case 'CYTBD'
                % 0.5o2 + q8h2 + 2he ---> q8 + h2o + 2he
                K4c = 0.449;
                K4 = 36; %h-1
                KUbiOxi = 6.5; %mmol/h.umol
                KmUbiOxi = 0.2; %in %
                EUbiOxi = 3;
%                 rUbiOxi_t = K4c*mu + K4*(1/(1+(Km4/pO2)^E3))*(1/(1+(pO2/Km4i)^E3i));
                UbiOxi_t = 0.5;
                UbiOxi_red = UbiOxi_t/(mc(q8)/mc(q8h2)*exp((pmf-EUbiOxi+Eq8)/Z)+1);
                flux(tfr) = KUbiOxi*UbiOxi_red*mc(q8h2)/(1+(KmUbiOxi/pO2)^EUbiOxi);
            case 'SUCDi'
                % succ + q8 <---> fum + q8h2
                flux(tfr) = Ksucdi*(Eq8-Z/2*log(mc(prid)/mc(sbid))-31);
            case 'FRD7'
                % succ + q8 <---> fum + q8h2
                flux(tfr) = Ksucdi*(Eq8-Z/2*log(mc(prid)/mc(sbid))-31);
        end
    end
end
                
            