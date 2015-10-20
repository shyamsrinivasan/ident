% Formulate the initial conecntration
function [mc,variable] = iconcentration(model,VMC,mc,assignFlag)
if nargin<4
    assignFlag = zeros(length(model.mets),1);
end
if nargin<3
    mc = zeros(model.nt_metab,1);
end
mets = fieldnames(VMC);
for jmc = 1:length(mets)
    tfm = strcmpi(model.mets,[mets{jmc} '[e]']);
    if any(tfm) && ~assignFlag(tfm)
        mc(tfm,:) = VMC.(mets{jmc});
        assignFlag(tfm,:) = 1;
    end
end
variable.MC = mc;