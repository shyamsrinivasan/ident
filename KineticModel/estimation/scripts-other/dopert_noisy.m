% generate perturbation data and add noise
% perturb system from noisy initial conditions
opts.x0 = noisy_xss(:,1);
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
pt_val = struct('exp_pid',{11,13,11,13},...
                'exp_pval',{[.5;1.0;1.5;2],...
                            [.1;.5;1.0;1.5;2],...
                            [.5;2],...
                            [.5;2]}); 
sol = getperturbations(pt_val,@perturb_nonoise,opts);

% add noise to perturbed data
pt_sol_id = 3;
noisy_sol = sol(pt_sol_id);
[noisy_sol.xss,noisy_sol.fss] = addnoise_wrapper(pt_sol_id,sol);