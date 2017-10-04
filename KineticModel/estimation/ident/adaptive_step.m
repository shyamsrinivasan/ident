% function for adpative step change in parameter for identifiability
% analysis
% obj - cas function and not a casadi symbolic object
% theta_k - value of optimized parameters from previous iteration
function [theta_step,obj_new,iter] =...
        adaptive_step(obj_k,theta_k,prob,p_val,fixed_pvalue,...
                      step_opts,input_data,ynoise_var)

% type = +1 for positive step and -1 for negative step
if isfield(step_opts,'type')
    type = step_opts.type;
else
    type = 1;
end
if isfield(step_opts,'PLE_threshold')
    delta_alpha = step_opts.PLE_threshold;
end
if isfield(step_opts,'minmax_step')
    min_step = step_opts.minmax_step(1);
    max_step = step_opts.minmax_step(2);
end
if isfield(prob,'xdynfun')
    xdynfun = prob.xdynfun; 
end
if isfield(prob,'xinit')
    xinit = prob.xinit; 
end
% if isfield(prob,'yinit')
%     yinit = prob.yinit;
% end
% if isfield(prob,'x')
%     p = prob.x; 
% end
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
% if isfield(data,'ynoise_var')
%     % measurement noise/error => y ~ N(0,sigma^2);
%     ynoise_var = data.ynoise_var;
% end

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
[x_newval,y_newval] =...
xdynfun(xinit,...
        repmat(theta_k,1,npts),...
        new_thetai,...
        repmat(p_val(6:9)',1,npts),...
        input_data);
% add initial value
% x_newval = [xinit x_newval];
% y_newval = [yinit y_newval];
% choose only points present in experimental data
y_model_newval = y_newval([1 3 4 5],:);

% create nlsqopt objective function
y_error = ((yexp-y_model_newval).^2)./ynoise_var;
obj_new = full(.5*sum(sum(y_error)));
obj_diff = obj_new-obj_k-q*delta_alpha;

while obj_diff>=eps && iter<=maxiter
    % cange thetai by theta_step
    new_thetai = fixed_pvalue+type*theta_step;
    
    [x_newval,y_newval] =...
    xdynfun(xinit,...
            repmat(theta_k,1,npts),...
            new_thetai,...
            repmat(p_val(6:9)',1,npts),...
            input_data);    
    % add initial value
%     x_newval = [xinit x_newval];
%     y_newval = [yinit y_newval];
    % choose only points present in experimental data
    y_model_newval = y_newval([1 3 4 5],:);

    % create nlsqopt objective function
    y_error = ((yexp-y_model_newval).^2)./ynoise_var;
    obj_new_sym = .5*sum(sum(y_error));
    obj_new = full(obj_new_sym);
    
    % calculate distance from threshold
    obj_diff = obj_new-obj_k-q*delta_alpha;
    
    theta_step = theta_step/1.1;
    iter = iter+1;
end

% lb < theta_step < ub to avoid PL discontinuities
if theta_step>max_step
    theta_step = max_step;
elseif theta_step<min_step
    theta_step=min_step;
end   