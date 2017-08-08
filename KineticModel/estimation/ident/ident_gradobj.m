function gradobj = ident_gradobj(x,p,data)

if isfield(data,'idx')
    idx = data.idx;
end
if isfield(data,'xexp')
    yexp = data.xexp;
end
if isfield(data,'nvar')
    nvar = data.nstate;
end
if isfield(data,'np')
    np = data.np;
end
if isifield(data,'xexp_var')
    yexp_var = data.yexp_var;
end
% function to solve ode and variational equations simultaneously using
% casadi
if isfield(data,'modelfh')
    gradfh = data.gradfh;
end

% augment parameter vector with fixed parameter vector value
aug_p = [x(1:idx-1);p(idx);x(idx:end)];

yres = gradfh(aug_p);
ysens = yres(nvar+1:end,:);
sens_mat = reshape(ysens,[np,nvar]);
sens_mat_copy = sens_mat;
sens_mat_copy(:,idx) = [];
