function [model,vl,vu] = changebounds(model,bounds,vl,vu,Vup_struct,fixgrowth)
if nargin < 6
    fixgrowth = 0;
end

%ETC
% vl(strcmpi(model.rxns,'ATPS4r')) = -100;
% vu(strcmpi(model.rxns,'ATPS4r')) = 0;
vl(strcmpi(model.rxns,'NADTRHD')) = 0;
vu(strcmpi(model.rxns,'NADTRHD')) = 100;
model.rev(strcmpi(model.rxns,'NADTRHD')) = 0;
vl(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'SUCCt2_2')) = 0;
vu(strcmpi(model.rxns,'FORt2')) = 0;

%Uptake Fluxes
if isfield(bounds,'Vuptake')
    if any(bounds.Vuptake)
        rxns = fieldnames(Vup_struct);
        for irxn = 1:length(rxns)
            tfr = strcmpi(model.rxns,rxns{irxn});
            if any(tfr)
                vl(tfr) = -bounds.Vuptake(tfr);
            end
        end
    end
%     if ~isempty(bounds.Vuptake)        
%         vl(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
%         vl(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
%     end
end

%Growth Fluxes
if fixgrowth
    if prxnid ~= model.bmrxn
        if isfield(model,'gmax')
            vl(model.bmrxn) = model.gmax;
            vu(model.bmrxn) = model.gmax;
        end
    end
else
    vl(model.bmrxn) = 0;    
end