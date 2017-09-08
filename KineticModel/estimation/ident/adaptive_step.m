% function for adpative step change in parameter for identifiability
% analysis
% obj - cas function and not a casadi symbolic object
% theta_k - value of optimized parameters from previous iteration
function [theta_step,obj_new,iter] =...
        adaptive_step(obj_k,theta_k,prob,p_val,fixed_pvalue,delta_alpha,freq,type)
if nargin<8
    type = 1;
end
% type = +1 for positive step and -1 for negative step

if isfield(prob,'xdynfun')
    xdynfun = prob.xdynfun; 
end
if isfield(prob,'xinit')
    xinit = prob.xinit; 
end
if isfield(prob,'x')
    p = prob.x; 
end
if isfield(prob,'npts')
    npts = prob.npts; 
end
if isfield(prob,'xexp') 
    xexp = prob.xexp; 
end
if isfield(prob,'yexp')
    yexp = prob.yexp;
else
    yexp = xexp;
end

maxiter = 1000;
% set threshold delta_alpha for theta_step
% conidition of new step in thetai
q = .1;  
iter = 1;
eps = 1e-4;

if type>0
    theta_step = 1;
elseif type<0
    theta_step = fixed_pvalue;
end

% first iteration
% cange thetai by theta_step
new_thetai = fixed_pvalue+type*theta_step;

% x_newval is not symbolic class(x_newval) = casadi.DM
x_newval =...
xdynfun(xinit,repmat(theta_k,1,npts),new_thetai,repmat(p_val(14:16)',1,npts),.1);
y_newval = x_newval;
% add initial value
y_newval = [xinit y_newval];
% choose only points present in experimental data
y_model_newval = y_newval(1:2,freq);

% create nlsqopt objective function
y_error = (yexp-y_model_newval);
obj_new = full(.5*dot(y_error,y_error));
obj_diff = obj_new-obj_k-q*delta_alpha;

while obj_diff>=eps && iter<=maxiter
    % cange thetai by theta_step
    new_thetai = fixed_pvalue+type*theta_step;
    
    x_newval =...
    xdynfun(xinit,repmat(theta_k,1,npts),new_thetai,repmat(p_val(14:16)',1,npts),.1);
    y_newval = x_newval;
    % add initial value
    y_newval = [xinit y_newval];
    % choose only points present in experimental data
    y_model_newval = y_newval(:,freq);

    % create nlsqopt objective function
    y_error = (xexp-y_model_newval);
    obj_new = full(.5*dot(y_error,y_error));
    
    % calculate distance from threshold
    obj_diff = obj_new-obj_k-q*delta_alpha;
    
    theta_step = theta_step/1.1;
    iter = iter+1;
end

% temp fixed step size
% theta_step = .01;




