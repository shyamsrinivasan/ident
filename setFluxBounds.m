function setFluxBounds(model,pvec)

%extracellular metabolites in M moles/L
met.glc = 0.2;
% met.h_e = 1e-7;
% met.h_c = 1e-7;
met.h2o_c = 55.0;%1.53e-13;
met.h2o_e = 50.0;%55.0;
met.o2_e = 0.0025;
met.pi_e = 1e-2;
% met.pi_c = 1e-3;
met.co2_e = 0.002;%1e-8;

[mc,assignFlag] = iconcentration(model,met);

nmetab = model.nint_metab;
next_metab = model.next_metab;
S = model.S;
K = pvec.K;

for im = 1:next_metab
    mid = im+nmetab;
    rxn = Vex(logical(S(mid,Vex)));
    
    if model.rev(rxn)
        
    elseif ~model.rev(rxn)
    end    
end