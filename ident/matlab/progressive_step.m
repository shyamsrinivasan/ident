function progressive_step(obj_last,theta_last,prob,p_val,fixed_pvalue,threshold)
if isfield(prob,'xdyn_fun')
    xdyn_fun = prob.xdyn_fun;
end
if isfield(prob,'objfun')
    objfun = prob.objfun;
end
if isfield(prob,'grad')
    grad = prob.grad;
end
if isfield(prob,'hess')
    hess = prob.hess;
end
if isfield(prob,'p')
    p = prob.p;
end
if isfield(prob,'ident_c')
    ident_c = prob.ident_c;
end
if isfield(prob,'p_useless')
    p_useless = prob.p_useless;
end
if isfield(prob,'acetate')
    acetate = prob.acetate;
end
if isfield(prob,'xexp') 
    xexp = prob.xexp; 
end
if isfield(prob,'xinit')
    xinit = prob.xinit;
end
if isfield(prob,'npts')
    npts = prob.npts; 
end

x_newval =...
xdynfun(xinit,repmat(p,1,npts),ident_c,repmat(p_useless,1,npts),acetate);
% add initial value
x_newval = [xinit x_newval];
% choose only points present in experimental data
x_model_newval = x_newval(:,1:100:2001);
% create nlsqopt objective function
x_error = (xexp-x_model_newval);

% theta_dash
theta_dash = 

% objective
obj_theta_last = full(.5*dot(x_error,x_error));
% obj_theta_last_fun = casadi.Function('obj_theta_last_fun',...
%                     {xinit,p,ident_c,p_useless,acetate},...
%                     {obj_theta_last});
% gradient
grad_theta_last = jacobian(obj_theta_last,p);
% grad_theta_last_fun = casadi.Function('grad_theta_last_fun',...
%                     {xinit,p,ident_c,p_useless,acetate},...
%                     {grad_theta_last});
% hessian                
hess_theta_last = hessian(obj_theta_last,p);
% hess_theta_last_fun = casadi.Function('hess_theta_last_fun',...
%                     {xinit,p,ident_c,p_useless,acetate},...
%                     {hess_theta_last});

chi2_theta_dash = obj_theta_last + grad_theta_last*

gamma = obj_theta_last_fun(xinit,theta_last,fixed_pvalue,p_val(14:16)',.1);
beta = hess_theta_last_fun(xinit,theta_last,fixed_pvalue,p_val(14:16)',.1);
theta_dash = theta - theta_last;

gamma = objfun();


obj_last = obj_theta_last +...
            gradient(theta_last)*theta_dash +...
            theta_dash'*hessian(theta_last)*theta_dash;

