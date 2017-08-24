function [exp_sol_start,noisy_sol_start] =...
        dopert_noisy(opts_temp,noisy_xss,noisy_fss,odep_bkp,pt_sol_id)
% generate perturbation data and add noise
% perturb system from noisy initial conditions

% create optionss structure for parallel loop
nstart = size(noisy_xss,2);
options = struct();
for ist = 1:nstart
    options(ist).x0 = opts_temp.x0;
    options(ist).tspan = opts_temp.tspan;
    options(ist).odep = opts_temp.odep;
    
    options(ist).x0 = noisy_xss(:,ist);
    options(ist).tspan = 0:.1:300;
    options(ist).odep = odep_bkp;    
end
% perturbations
pt_val = struct('exp_pid',{11,12,13},...
                'exp_pval',{[2],...
                            [.2],...
                            [.5]}); 
% id = 1;
noisy_sol_start = cell(nstart,1);
exp_sol_start = cell(nstart,1);
if nstart>10
    parfor ist = 1:nstart 
        xss_temp = noisy_xss(:,ist);
        fss_temp = noisy_fss(:,ist);

        opts_temp = options(ist);
        sol = getperturbations(pt_val,@perturb_nonoise,opts_temp);
        % close all

        % add noise to perturbed data
        % pt_sol_id = 3;
        noisy_sol = sol(pt_sol_id);
        noisy_ss = addnoise_wrapper(pt_sol_id,sol);
        for j = 1:length(pt_sol_id)
            noisy_sol(j).xss = noisy_ss(j).xss;
            noisy_sol(j).fss = noisy_ss(j).fss;
        end
        % include wt (unperturbed) data
        noisy_sol(j+1).xss = xss_temp;
        noisy_sol(j+1).fss = fss_temp;
        noisy_sol(j+1).exp_pid = 0;
        noisy_sol(j+1).exp_pval = 0;
        noisy_sol(j+1).odep = odep_bkp;

        % combine all perturbation flux data into a single vector
        exp_sol_start{ist}.xss = cat(2,noisy_sol.xss);
        exp_sol_start{ist}.fss = cat(2,noisy_sol.fss);
        exp_sol_start{ist}.exp_pid = cat(2,noisy_sol.exp_pid);
        exp_sol_start{ist}.exp_pval = cat(2,noisy_sol.exp_pval);
        exp_sol_start{ist}.odep = cat(1,noisy_sol.odep);

        noisy_sol_start{ist} = noisy_sol;
    %     exp_sol_start{ist} = exp_sol;
    end
else
    for ist = 1:nstart 
        xss_temp = noisy_xss(:,ist);
        fss_temp = noisy_fss(:,ist);

        opts_temp = options(ist);
        sol = getperturbations(pt_val,@perturb_nonoise,opts_temp);
        % close all

        % add noise to perturbed data
        % pt_sol_id = 3;
        noisy_sol = sol(pt_sol_id);
        noisy_ss = addnoise_wrapper(pt_sol_id,sol);
        for j = 1:length(pt_sol_id)
            noisy_sol(j).xss = noisy_ss(j).xss;
            noisy_sol(j).fss = noisy_ss(j).fss;
        end
        % include wt (unperturbed) data
        noisy_sol(j+1).xss = xss_temp;
        noisy_sol(j+1).fss = fss_temp;
        noisy_sol(j+1).exp_pid = 0;
        noisy_sol(j+1).exp_pval = 0;
        noisy_sol(j+1).odep = odep_bkp;

        % combine all perturbation flux data into a single vector
        exp_sol_start{ist}.xss = cat(2,noisy_sol.xss);
        exp_sol_start{ist}.fss = cat(2,noisy_sol.fss);
        exp_sol_start{ist}.exp_pid = cat(2,noisy_sol.exp_pid);
        exp_sol_start{ist}.exp_pval = cat(2,noisy_sol.exp_pval);
        exp_sol_start{ist}.odep = cat(1,noisy_sol.odep);

        noisy_sol_start{ist} = noisy_sol;
    %     exp_sol_start{ist} = exp_sol;
    end    
end


