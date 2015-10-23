% Formulate the initial conecntration
function [mc,assignFlag] = iconcentration(model,VMC,mc,assignFlag)
if nargin<4
    assignFlag = zeros(length(model.mets),1);
end
if nargin<3
    mc = zeros(model.nt_metab,1);
end

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
    if any(tfm) %&& ~assignFlag(tfm)
        mc(tfm,:) = VMC.(mets{jmc});
        assignFlag(tfm,:) = 1;
    end
end
