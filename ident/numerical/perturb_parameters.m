% perform parameter perturbations based on info in perturbation srtcut
% array
function solution = perturb_parameters(odefun,solver_opts,...
                                       original_parameters,...
                                       perturbation,flux_fun)
if nargin<5
    flux_fun = [];
end
number_perturbations = size(perturbation, 2);

solution = struct([]);
for i_pert = 1:number_perturbations
    select_parameters = original_parameters;
    parameter_id = perturbation(i_pert).id;
    parameter_change = perturbation(i_pert).change;
    % change parameter value
    select_parameters(parameter_id) =...
        select_parameters(parameter_id)*(1+parameter_change);
    one_solution = run_ode(odefun, solver_opts, select_parameters, flux_fun); 
    solution(i_pert).t = one_solution.t;
    solution(i_pert).y = one_solution.y;
    solution(i_pert).yss = one_solution.yss;
    if isfield(one_solution, 'flux')
        solution(i_pert).flux = one_solution.flux;
    else
        solution(i_pert).flux = [];
    end
    if isfield(one_solution, 'fss')
        solution(i_pert).fss = one_solution.fss;
    else
        solution(i_pert).fss = [];
    end
    solution(i_pert).parameter = select_parameters;
end