function noisy_sol = dopert_noisy(opts,noisy_xss,odep_bkp,pt_sol_id)
% generate perturbation data and add noise
% perturb system from noisy initial conditions
opts.x0 = noisy_xss(:,1);
opts.tspan = 0:.1:300;
opts.odep = odep_bkp;
pt_val = struct('exp_pid',{11,12,13},...
                'exp_pval',{[.5;1;2],...
                            [.5;1;2],...
                            [.5;1;2]}); 
sol = getperturbations(pt_val,@perturb_nonoise,opts);
% close all

% add noise to perturbed data
% pt_sol_id = 3;
noisy_sol = sol(pt_sol_id);
noisy_ss = addnoise_wrapper(pt_sol_id,sol);
for j = 1:length(pt_sol_id)
    noisy_sol(j).xss = noisy_ss(j).xss;
    noisy_sol(j).fss = noisy_ss(j).fss;
end
