% run the looping construct of the PLE algorithm for identifiability
function [PLEvals] =...
getPLE(thetai_fixed_value,theta_step,p0,p_val,delta_alpha,data,maxiter,pos_neg)
if nargin<8
    pos_neg = 2;
end
if nargin<7
    maxiter = 100;
end

if isfield(data,'plim')
    plim = data.plim;
end

% positive step
chiPLE_pos = zeros(1,maxiter);
xPLE_pos = zeros(length(p0),maxiter);
thetai_inc_pos = zeros(1,maxiter);
thetai_step_pos = zeros(1,maxiter);
obj_step_pos = zeros(1,maxiter);

thetai_fixed_value_start = thetai_fixed_value;
theta_step_start = theta_step;
p0_start = p0;

if pos_neg==2 || pos_neg==1
    iter_pos = 1;
    pre_chiPLE = chiPLE_pos(1);
    while iter_pos<=maxiter && pre_chiPLE<=delta_alpha &&...
                            thetai_fixed_value<=plim(2)    

        [optsol,thetai_fixed_value,theta_step,obj_new] =...
        PLEiter(thetai_fixed_value,theta_step,p0,p_val,delta_alpha,data,1);

        % new initial p0 = old optimal value
        p0 = optsol.xval;    

        % store obj values
        obj_step_pos(iter_pos) = obj_new;
        chiPLE_pos(iter_pos) = optsol.fval;
        xPLE_pos(:,iter_pos) = optsol.xval;
        thetai_inc_pos(iter_pos) = thetai_fixed_value;
        thetai_step_pos(iter_pos) = theta_step;

        pre_chiPLE = optsol.fval;    
        iter_pos = iter_pos+1;
    end
else
    iter_pos = 1;
end

% negative step
chiPLE_neg = zeros(1,maxiter);
xPLE_neg = zeros(length(p0),maxiter);
thetai_inc_neg = zeros(1,maxiter);
thetai_step_neg = zeros(1,maxiter);
obj_step_neg = zeros(1,maxiter);

if pos_neg==2 || pos_neg==3
    thetai_fixed_value = thetai_fixed_value_start;
    theta_step = theta_step_start;
    p0 = p0_start;
    iter_neg = 1;
    pre_chiPLE = chiPLE_neg(1);
    while iter_neg<=maxiter && pre_chiPLE<=delta_alpha &&...
                            thetai_fixed_value>=plim(1)    

        [optsol,thetai_fixed_value,theta_step,obj_new] =...
        PLEiter(thetai_fixed_value,theta_step,p0,p_val,delta_alpha,data,-1);

        % new initial p0 = old optimal value
        p0 = optsol.xval;    

        % store obj values
        obj_step_neg(iter_neg) = obj_new;
        chiPLE_neg(iter_neg) = optsol.fval;
        xPLE_neg(:,iter_neg) = optsol.xval;
        thetai_inc_neg(iter_neg) = thetai_fixed_value;
        thetai_step_neg(iter_neg) = theta_step;

        pre_chiPLE = optsol.fval;    
        iter_neg = iter_neg+1;
    end
else
    iter_neg = 1;
end

% calculate CI based on results from above

% collate results of both pos and neg steps on thetai
PLEvals = struct();
PLEvals.obj_step = [fliplr(obj_step_neg(1:iter_neg-1)) obj_step_pos(1:iter_pos-1)];
PLEvals.chiPLE = [fliplr(chiPLE_neg(1:iter_neg-1)) chiPLE_pos(1:iter_pos-1)];
PLEvals.xPLE = [fliplr(xPLE_neg(:,1:iter_neg-1)) xPLE_pos(:,1:iter_pos-1)];
PLEvals.thetai_inc = [fliplr(thetai_inc_neg(1:iter_neg-1)),... 
                        thetai_inc_pos(1:iter_pos-1)];
PLEvals.thetai_step = [fliplr(thetai_step_neg(1:iter_neg-1)),...
                        thetai_step_pos(1:iter_pos-1)];
PLEvals.iter = [iter_pos-1 iter_neg-1];

