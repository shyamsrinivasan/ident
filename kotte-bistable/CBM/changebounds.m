% function [model,vl,vu] = changebounds(model,bounds,vl,vu,Vup_struct,fixgrowth)
function [model,bounds] = changebounds(model,ess_rxn,bounds,fixgrowth)
if nargin < 4
    fixgrowth = 0;
end
model.rev(strcmpi(model.rxns,'NADTRHD')) = 0;
if isfield(model,'Vuptake')
    Vuptake = model.Vuptake;
else
    Vuptake = [];
end
if nargin < 3
    bounds = struct();    
    nr = size(model.S,2);
    vl = zeros(nr,1);        
    vl(logical(model.rev)) = -100;
    vu = zeros(nr,1);          
    vu(vu==0) = 100;
else    
    vl = bounds.vl;
    vu = bounds.vu;    
end
if nargin < 2
    ess_rxn = {};
end
% model.rev(strcmpi(model.rxns,'NADTRHD')) = 0;

% Vuptake = model.Vuptake;
% vl = zeros(model.nt_rxn,1);        
% vl(logical(model.rev)) = -100;
% vu = zeros(model.nt_rxn,1);          
% vu(vu==0) = 100;

%set other flux bounds
%ETC
% vl(strcmpi(model.rxns,'ATPS4r')) = -100;
% vu(strcmpi(model.rxns,'ATPS4r')) = 0;
vl(strcmpi(model.rxns,'NADTRHD')) = 0;
vu(strcmpi(model.rxns,'NADTRHD')) = 500;
vl(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'THD2')) = 0;
vu(strcmpi(model.rxns,'SUCCt2_2')) = 0;
vu(strcmpi(model.rxns,'FORt2')) = 0;
vl(strcmpi(model.rxns,'PPC')) = 0;

vu(strcmpi(model.rxns,'ATPS4r')) = 200;


%Uptake Fluxes
if ~isempty(Vuptake) && any(Vuptake)
    vl(logical(Vuptake)) = -Vuptake(logical(Vuptake));
%     vu(logical(Vuptake)) = -Vuptake(logical(Vuptake));
end

%change bounds for exchange metabolites
% ess_rxn = {'exCO2','exH','exH2O','exPI','exO2','exGLC'};
essid = [];
for iess = 1:length(ess_rxn)
    essid = union(essid,find(strcmpi(ess_rxn{iess},model.rxns)));
end
if isfield(model,'VFex')
    Vess = setdiff(model.VFex,essid);
    vl(Vess) = 0;
end

%atp maintanance
vl(strcmpi(model.rxns,'ATPM')) = 8.39;
vu(strcmpi(model.rxns,'ATPM')) = 200;

%Growth Fluxes
if fixgrowth
    if prxnid ~= model.bmrxn
        if isfield(model,'gmax')
            vl(model.bmrxn) = model.gmax;
            vu(model.bmrxn) = model.gmax;
        else
            error('Nogmax:changebounds','No field gmax found in model. Cannot fix growth');
        end
    end
else
    vl(model.bmrxn) = 0;    
end

bounds.Vuptake = Vuptake;
bounds.vl = vl;
bounds.vu = vu;

% %Growth Fluxes - called separately to be set
% if nargin == 3
%     if vl
%         if isfield(model,'gmax')
%             bounds.vl(model.bmrxn) = model.gmax;
%             bounds.vu(model.bmrxn) = model.gmax;
%         end
%     else
%         bounds.vl(model.bmrxn) = 0;
%     end
%     vl = bounds.vl;
%     vu = bounds.vu;
%     return
% end
% 
% if nargin < 6
%     fixgrowth = 0;
% end
% %ETC
% % vl(strcmpi(model.rxns,'ATPS4r')) = -100;
% % vu(strcmpi(model.rxns,'ATPS4r')) = 0;
% vl(strcmpi(model.rxns,'NADTRHD')) = 0;
% vu(strcmpi(model.rxns,'NADTRHD')) = 100;
% model.rev(strcmpi(model.rxns,'NADTRHD')) = 0;
% vl(strcmpi(model.rxns,'THD2')) = 0;
% vu(strcmpi(model.rxns,'THD2')) = 0;
% vu(strcmpi(model.rxns,'SUCCt2_2')) = 0;
% vu(strcmpi(model.rxns,'FORt2')) = 0;

%Uptake Fluxes
% if isfield(bounds,'Vuptake')
%     if any(bounds.Vuptake)
%         vl(logical(bounds.Vuptake)) = bounds.Vuptake(logical(bounds.Vuptake));
% %         rxns = fieldnames(Vup_struct);
% %         for irxn = 1:length(rxns)
% %             tfr = strcmpi(model.rxns,rxns{irxn});
% %             if any(tfr)
% %                 vl(tfr) = -bounds.Vuptake(tfr);
% %             end
% %         end
%     end
% %     if ~isempty(bounds.Vuptake)        
% %         vl(strcmpi(model.rxns,'exGLC')) = -bounds.Vuptake(strcmpi(model.rxns,'exGLC'));        
% %         vl(strcmpi(model.rxns,'exO2')) = -bounds.Vuptake(strcmpi(model.rxns,'exO2'));
% %     end
% end

%Growth Fluxes
% if fixgrowth
%     if prxnid ~= model.bmrxn
%         if isfield(model,'gmax')
%             vl(model.bmrxn) = model.gmax;
%             vu(model.bmrxn) = model.gmax;
%         end
%     end
% else
%     vl(model.bmrxn) = 0;    
% end