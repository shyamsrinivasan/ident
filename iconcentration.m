% Formulate the initial conecntration
function [varargout] = iconcentration(model,VMC,mc_lb,assignFlag)
if nargin<4
    assignFlag = zeros(length(model.mets),1);
end
if nargin<3
    mc_lb = zeros(model.nt_metab,1);
    mc_ub = zeros(model.nt_metab,1);
else
    mc_ub = mc_lb;
end
% mc_ub = zeros(model.nt_metab,1);

if ~isempty(VMC)
    mets = fieldnames(VMC);
    kmet = cell(length(mets),1);
    for jmc = 1:length(mets)
        ke = strfind(mets{jmc},'_e');
        kc = strfind(mets{jmc},'_c');
        if ~isempty(ke)
            kmet{jmc} = mets{jmc}(1:ke-1);
            cmp = '[e]';
        elseif ~isempty(kc)
            kmet{jmc} = mets{jmc}(1:kc-1);
            cmp = '[c]';
        end
        
        tfm = strcmpi(model.mets,[kmet{jmc} cmp]);    
        if any(tfm) && ~assignFlag(tfm)
            if length(VMC)>1
                mc_lb(tfm,:) = VMC(1).(mets{jmc});
                mc_ub(tfm,:) = VMC(2).(mets{jmc});
                assignFlag(tfm,:) = 1;
            else
                mc_lb(tfm,:) = VMC.(mets{jmc});
                assignFlag(tfm,:) = 1;
            end
        end
    end
end
if ~any(mc_ub)
    mc_ub = [];
end
varargout{1} = mc_lb;
varargout{2} = mc_ub;
if mc_ub == mc_lb
    varargout{2} = logical(assignFlag);
else
    varargout{3} = logical(assignFlag);
end
