% objective for identifiability analysis
% sum of square minimization of variance weighted error
function objfx = ident_obj(x,data)

if isfield(data,'odep')
    p = data.odep;
end
if isfield(data,'idx')
    idx = data.idx;
end
if isfield(data,'xexp_dyn')
    yexp = data.xexp_dyn;
end
if isfield(data,'xexp_var')
    yexp_var = data.xexp_var;
end
if isfield(data,'modelf')
    modelf = data.modelf;
end

% augment parameter vector with fixed parameter vector value
aug_p = [x(1:idx-1) p(idx) x(idx:end)];

% calculate ymodel(theta,t) with augmented parameter vector
ymodel = modelf(aug_p);

objfx = sum(sum(((yexp-ymodel)./yexp_var).^2,2));


                              






