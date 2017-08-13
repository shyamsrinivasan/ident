% function for adpative step change in parameter for identifiability
% analysis
% obj - cas function and not a casadi symbolic object
% theta_k - value of optimized parameters from previous iteration
function adaptive_step(obj_k,theta_k,prob,p_val,fixed_pvalue)

if isfield(prob,'xdynfun'), xdynfun = prob.xdynfun; end
if isfield(prob,'xinit'), xinit = prob.xinit; end
if isfield(prob,'x'), p = prob.x; end
if isfield(prob,'npts'), npts = prob.npts; end
if isfield(prob,'xexp'), xexp = prob.xexp; end

theta_step = .001;


new_thetai = fixed_pvalue+theta_step;
x_newval =...
xdynfun(xinit,repmat(theta_k,1,npts),new_thetai,repmat(p_val(14:16)',1,npts),.1);
% add initial value
x_newval = [xinit x_newval];
% choose only points present in experimental data
x_model_newval = x_newval(:,1:100:2001);

% create nlsqopt objective function
x_error = (xexp-x_model_newval);
obj = .5*dot(x_error,x_error);

% conidition of new step in thetai
obj_k-obj-q*delta_alpha <= eps;




