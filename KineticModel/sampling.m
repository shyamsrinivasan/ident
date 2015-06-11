%Sample Quantities to build ensemble of models
function [ensb] = build_ensemble(nmodels,model)
nmodels = 100;
nsamples = nmodels;
nrxn = model.nrxns;

for isample = 1:nsamples
    mname = sprintf('model%d',isample);
    ensb.(mname).norm_metab = model.S;
    ensb.(mname).norm_reg = model.SI;
end   

metab = cell(nrxn,1);
nmetab = zeros(nrxn,1);
reg = cell(nrxn,1);
nreg = zeros(nrxn,1);
for irxn = 1:nrxn
    metab{irxn} = find(model.S(:,irxn));
    nmetab(irxn) = length(metab{irxn});
    reg{irxn} = find(model.SI(:,irxn));
    nreg(irxn) = length(reg{irxn});
end
        
for isample = 1:nsamples
    mname = sprintf('model%d',isample);    
    for irxn = 1:nrxn        
        %Substrate/Product Samples            
        met_sat = betarnd(1.5,4.5,nmetab(irxn),1); %Beta Distribution 
        ensb.(mname).norm_metab(metab{irxn}) = met_sat./(1-met_sat);
        %Regulator Samples             
        reg_sat = betarnd(1.5,4.5,nreg(irxn),1);
        ensb.(mname).norm_reg(reg{irxn}) = reg_sat./(1-reg_sat);
    end   
end