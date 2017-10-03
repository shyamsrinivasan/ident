function [optsol,thetai_fixed_value,theta_step,obj_new,iter_theta_step] =...
        PLEiter(thetai_fixed_value,theta_step,xval,p_val,PLE_threshold,setup_opts,type)
    
% new theta
thetai_fixed_value = thetai_fixed_value + type*theta_step;

% PLE through optimization for new thetai
% prob_cas = identopt_setup(setup_opts,thetai_fixed_value);
prob_cas = setpleopt_ss(setup_opts,thetai_fixed_value);
optsol = solve_nlsqopt(prob_cas,xval);

% adpative step in thetai
adap_step_opts.PLE_threshold = PLE_threshold;
adap_step_opts.minmax_step = setup_opts.minmax_step;
adap_step_opts.type = type;
[theta_step,obj_new,iter_theta_step] =...
adaptive_step(optsol.fval,optsol.xval,prob_cas,p_val,...
              thetai_fixed_value,adap_step_opts,setup_opts.freq);
    
    