function gradobj = ident_gradobj(x,data,funh)

if isfield(funh,'modelsensfh')
    modelsensfh = funh.modelsensfh;
end
if isfield(data,'odep')
    p = data.odep;
end
if isfield(data,'idx')
    idx = data.idx;
end
if isfield(data,'nvar')
    nvar = data.nstate;
end
if isfield(data,'np')
    np = data.np;
end
% if isfield(data,'xexp')
%     yexp = data.xexp;
% end
% if isifield(data,'xexp_var')
%     yexp_var = data.yexp_var;
% end
% if isfield(data,'gradobj')
%     gradobjfh = data.gradobj;
% end


% augment parameter vector with fixed parameter vector value
aug_p = [x(1:idx-1);p(idx);x(idx:end)];

% calculate sensitivity by solving augmented system
yres = modelsensfh(aug_p);
ymodel = yres(1:nvar,:);
ysens = yres(nvar+1:end,:);
sens_mat = reshape(ysens,[nvar,np]);
sens_mat_copy = sens_mat;
sens_mat_copy(:,idx) = [];




% function to solve ode and variational equations simultaneously using
% casadi
if isfield(data,'modelfh')
    gradfh = data.gradfh;
end