% test if x0 satisfies constraints 
function [eqcons,lecons,gecons] = checkconstr(consfh,xval,nlrhs,nle)

nlconsval = consfh(xval);
% == constraints
eqcons = ones(length(nlrhs),1);
if ~all(nlconsval(nle==0) == nlrhs(nle==0))
    eqcons(nlconsval ~= nlrhs) = -1;
    eqcons(nle~=0) = 0;    
end
% <= constraints
lecons = ones(length(nlrhs),1);
if ~all(nlconsval(nle==-1) == nlrhs(nle==-1))
    lecons(nlconsval > nlrhs) = -1;
    lecons(nle~=-1) = 0;    
end
% >= constraints
gecons = ones(length(nlrhs),1);
if ~all(nlconsval(nle==1) == nlrhs(nle==1))
    gecons(nlconsval < nlrhs) = -1;
    gecons(nle~=1) = 0;    
end