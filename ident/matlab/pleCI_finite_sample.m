% get finite sample confidence intervals (CI) for PLE 
% does not need hessians of the objective function at the optimal solution
function sigma =...
        pleCI_finite_sample(delta_alpha,chiPLE,thetai,thetai_start)

% theta s.t. obj(theta)-obj(optimal_theta)<delta_alpha
% get point of intersection of PL curve with threshold line (delta_alpha)
sigma = zeros(1,2);
pl_diff = delta_alpha-chiPLE;

% find position of sign change in pl_diff
bounds_lb_id = find(sign(pl_diff)<0,1,'first');
bounds_ub_id = find(sign(pl_diff)<0,1,'last');

if ~isempty(bounds_lb_id)
    bounds_lb = thetai(bounds_lb_id);
else
    bounds_lb = [];
end
if ~isempty(bounds_ub_id)
    bounds_ub = thetai(bounds_ub_id);
else
    bounds_ub = [];
end

if ~isempty(bounds_lb) && ~isempty(bounds_ub)
    if bounds_lb==bounds_ub % only one bound (practical non-identifiability)
        if bounds_ub>thetai_start % bound is ub
            sigma(2) = bounds_ub;
            sigma(1) = -Inf;
        elseif bounds_lb<thetai_start % bound is lb
            sigma(1) = bounds_lb;
            sigma(2) = +Inf;
        end
    else
        sigma(1) = bounds_lb;
        sigma(2) = bounds_ub;
    end
end



        









