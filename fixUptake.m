function Vuptake = fixUptake(model,Vup_struct)
%fix uptake of reactions designated as fields in Vuptake to model.Vuptake

if ~isfield(model,'Vuptake')
    Vuptake = zeros(model.nt_rxn,1);
    if ~isempty(Vup_struct)
        rxns = fieldnames(Vup_struct);
    else
        rxns = {};
    end
    if ~isempty(rxns)
        for irxn = 1:length(rxns)
            tfr = strcmpi(model.rxns,rxns{irxn});
            if any(tfr)
                Vuptake(tfr) = Vup_struct.(rxns{irxn});
            end
        end
    end
%     model.Vuptake = Vuptake;
end