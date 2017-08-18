% get finite sample confidence intervals (CI) for PLE 
% does not need hessians of the objective function at the optimal solution
function [sigma,collect_eps] = pleCI_finite_sample(delta_alpha,chiPLE,thetai)

% theta s.t. obj(theta)-obj(optimal_theta)<delta_alpha
% get point of intersection of PL curve with threshold line (delta_alpha)
eps = 1e-1;
pl_diff = delta_alpha-chiPLE;

while any(pl_diff<eps)    
    if ~any(pl_diff<eps)
        collect_eps = eps;
        collect_pos_first = find(pl_diff<eps,1,'first');
        collect_pos_last = find(pl_diff<eps,1,'last');
    else
        collect_eps = 0;
        collect_pos_first = false;
        collect_pos_last = false;
    end
    eps = eps/5;
end

sigma = zeros(1,2);
if collect_pos_first==collect_pos_last
    % only one confidence interval is finite
    sigma(1) = thetai(collect_pos_first);
    sigma(2) = sigma(1); %
else
    % finite confidence intervals on both ends
    % CI lb
    sigma(1) = thetai(collect_pos_first);
    sigma(2) = thetai(collect_pos_last);
end
        









