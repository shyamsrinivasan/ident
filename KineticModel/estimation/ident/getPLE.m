% run the looping construct of the PLE algorithm for identifiability
function [PLEvals] =...
getPLE(thetai_fixed_value,theta_step,p0,p_val,delta_alpha,data,maxiter)
if nargin<7
    maxiter = 100;
end

chiPLE = zeros(1,maxiter);
xPLE = zeros(length(p0),maxiter);
thetai_inc = zeros(1,maxiter);
thetai_step = zeros(1,maxiter);
obj_step = zeros(1,maxiter);

iter = 1;
pre_chiPLE = chiPLE(1);
while iter<=maxiter && pre_chiPLE<=delta_alpha
    
    [optsol,thetai_fixed_value,theta_step,obj_new] =...
    PLEiter(thetai_fixed_value,theta_step,p0,p_val,delta_alpha,data);
        
    % new initial p0 = old optimal value
    p0 = optsol.xval;    
     
    % store obj values
    obj_step(iter) = obj_new;
    chiPLE(iter) = optsol.fval;
    xPLE(:,iter) = optsol.xval;
    thetai_inc(iter) = thetai_fixed_value;
    thetai_step(iter) = theta_step;
    
    pre_chiPLE = optsol.fval;
    
    % figure for PLE
%     line(thetai_inc,chiPLE(iter),'LineStyle','none','Marker','.','MarkerSize',10);
    
    iter = iter+1;
end

PLEvals = struct();
PLEvals.obj_step = obj_step;
PLEvals.chiPLE = chiPLE;
PLEvals.xPLE = xPLE;
PLEvals.thetai_inc = thetai_inc;
PLEvals.thetai_step = thetai_step;


