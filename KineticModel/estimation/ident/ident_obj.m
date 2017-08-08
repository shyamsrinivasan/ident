% objective for identifiability analysis
% sum of square minimization of variance weighted error
function objfx = ident_obj(x,p,data)

if isfield(data,'idx')
    idx = data.idx;
end
if isfield(data,'xexp')
    yexp = data.xexp;
end
if isifield(data,'xexp_var')
    yexp_var = data.yexp_var;
end
if isfield(data,'modelfh')
    modelfh = data.modelfh;
end

% augment parameter vector with fixed parameter vector value
aug_p = [x(1:idx-1);p(idx);x(idx:end)];

% calculate ymodel(theta,t) with augmented parameter vector
ymodel = modelfh(aug_p);

objfx = sum((yexp-ymodel)./yexp_var).^2;


                              






