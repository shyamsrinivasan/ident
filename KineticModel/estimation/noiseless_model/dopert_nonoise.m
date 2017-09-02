function [exp_sol,no_noise_sol] =...
        dopert_nonoise(opts,xss,fss,odep_bkp,pt_sol_id)
    
% generate perturbation data 
% perturb system from deterministic (non-noisy) initial conditions
id = 1;
opts.x0 = xss(:,id);
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
pt_val = struct('exp_pid',{11,12,13,11,12,13},...
                'exp_pval',{[2],...
                            [.1],...
                            [.1],...
                            [.1],...
                            [2],...
                            [2]}); 
sol = getperturbations(pt_val,@perturb_nonoise,opts);

no_noise_sol = sol(pt_sol_id);
j = length(pt_sol_id);
% include wt (unperturbed) data
no_noise_sol(j+1).xss = xss(:,id);
no_noise_sol(j+1).fss = fss(:,id);
no_noise_sol(j+1).exp_pid = 0;
no_noise_sol(j+1).exp_pval = 0;
no_noise_sol(j+1).odep = odep_bkp;

% combine all perturbation flux data into a single vector
exp_sol.xss = cat(2,no_noise_sol.xss);
exp_sol.fss = cat(2,no_noise_sol.fss);
exp_sol.exp_pid = cat(2,no_noise_sol.exp_pid);
exp_sol.exp_pval = cat(2,no_noise_sol.exp_pval);
exp_sol.odep = cat(1,no_noise_sol.odep);

    
    
