% Formulate the initial conecntration
function [mc,variable] = iconcentration(model,VMC,mc)
if nargin<3
    mc = zeros(model.nt_metab,1);
end
mets = fiedlnames(VMC);
for jmc = 1:length(mets)
    tfm = strcmpi(model.mets,[mets{jmc} '[e]']);
    if any(tfm)
        mc(tfm) = VMC.(mets{jmc});
    end
end
variable.MC = mc;