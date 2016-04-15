function varargout = assignFlux(model,VMC,fx_lb,assignFlag)
if nargin<4
    assignFlag = zeros(length(model.rxns),1);
end
if nargin<3
    fx_lb = zeros(model.nt_rxn,1);
    fx_ub = zeros(model.nt_rxn,1);
else
    fx_ub = fx_lb;
end

if ~isempty(VMC)
    rxns = fieldnames(VMC);    
    for jrx = 1:length(rxns)
        tfr = strcmpi(model.rxns,rxns{jrx});
        if any(tfr)
            if length(VMC)>1
                fx_lb(tfr,:) = VMC(1).(rxns{jrx});
                fx_ub(tfr,:) = VMC(2).(rxns{jrx});
                assignFlag(tfr,:) = 1;
            else
                fx_lb(tfr,:) = VMC.(rxns{jrx});
                assignFlag(tfr,:) = 1;
            end
        end
    end
end
if ~any(fx_ub)
    fx_ub = [];
end

varargout{1} = fx_lb;
varargout{2} = fx_ub;
varargout{3} = assignFlag;