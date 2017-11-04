% function to test whether given perturbation will satisfy the conditiond
% for identifiability
% fh - function for testing identifiability of a specific flux
% np - # perturbations required for testing
function check_pert_ident(exp_sol,np,fh)

nexpt = size(exp_sol,2);
if nexpt>np
    % check all np combinations of perturbations in exp_sol
    for i = 1:nexpt
        for j = 1:nexpt
            if j~=i
                res = fh();
            end
        end
    end
elseif nexpt==np
    % pass all existing perturbations for testing
end
