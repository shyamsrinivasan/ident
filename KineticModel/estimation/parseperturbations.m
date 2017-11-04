function [exp_sol,newsol] = parseperturbations(sol,use_pert)

npert = size(sol,2);
npert_used = length(use_pert);
if npert_used<npert
    newsol = sol(use_pert);
    exp_sol.xss = cat(2,newsol.xss);
    exp_sol.fss = cat(2,newsol.fss);
    if isfield(newsol,'noise_xss')
        exp_sol.noise_xss = cat(2,newsol.noise_xss);
    end
    if isfield(newsol,'noise_fss')
        exp_sol.noise_fss = cat(2,newsol.noise_fss);
    end
    p_pert = cat(1,newsol.odep);
    p_wt = repmat(p_pert(end,11:13),npert_used,1);    
    p_pert = p_pert(:,11:13);    
    p_pert_logical = (p_pert~=p_wt);
    % convert all wt parameters to 1
    p_pert_logical(all(p_pert_logical(:,:)'==0)',:) = 1;
    exp_sol.p_pert_logical = p_pert_logical';
    exp_sol.p_pert = p_pert';
    p_pert = mat2cell(p_pert,ones(1,npert_used));
    p_pert_logical = mat2cell(p_pert_logical,ones(1,npert_used));
    [newsol(:).p_pert] = p_pert{:};
    [newsol(:).p_pert_logical] = p_pert_logical{:};
    if isfield(newsol,'xdyn')
        if ismember(npert,use_pert)
            % for wt - use final ss value all along the time horizon
            wt_xdyn = newsol(use_pert==npert).xdyn;
            wt_xdyn_ss = wt_xdyn(:,end);
            t_freq = size(wt_xdyn,2);
            newsol(use_pert==npert).xdyn = repmat(wt_xdyn_ss,1,t_freq);
        end        
        exp_sol.xdyn = cat(1,newsol.xdyn);
        if isfield(newsol,'noise_xdyn')
            exp_sol.noise_xdyn = cat(1,newsol.noise_xdyn);
        end
    end
    if isfield(newsol,'fdyn')
        exp_sol.fdyn = cat(1,newsol.fdyn);
        if isfield(newsol,'noise_fdyn')
            exp_sol.noise_fdyn = cat(1,newsol.noise_fdyn);
        end
    end
    exp_sol.exp_pid = cat(2,newsol.exp_pid);
    exp_sol.exp_pval = cat(2,newsol.exp_pval);
    exp_sol.odep = cat(1,newsol.odep);
else
    newsol = sol;
end