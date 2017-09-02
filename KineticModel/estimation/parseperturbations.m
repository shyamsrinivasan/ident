function [exp_sol,newsol] = parseperturbations(sol,use_pert)

if length(use_pert)<size(sol,2)
    newsol = sol(use_pert);
    exp_sol.xss = cat(2,newsol.xss);
    exp_sol.fss = cat(2,newsol.fss);
    if isfield(newsol,'xdyn')
        if ismember(size(sol,2),use_pert)
            % for wt - use final ss value all along the time horizon
            wt_xdyn = newsol(use_pert==size(sol,2)).xdyn;
            wt_xdyn_ss = wt_xdyn(:,end);
            t_freq = size(wt_xdyn,2);
            newsol(use_pert==size(sol,2)).xdyn = repmat(wt_xdyn_ss,1,t_freq);
        end        
        exp_sol.xdyn = cat(1,newsol.xdyn);
    end
    if isfield(newsol,'fdyn')
        exp_sol.fdyn = cat(1,newsol.fdyn);
    end
    exp_sol.exp_pid = cat(2,newsol.exp_pid);
    exp_sol.exp_pval = cat(2,newsol.exp_pval);
    exp_sol.odep = cat(1,newsol.odep);
else
    newsol = sol;
end