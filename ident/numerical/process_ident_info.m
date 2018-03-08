% process numerical ientifiability information after solving systems of
% nlae to identify parameters
function data = process_ident_info(ident_solution)

% collect only non-zero fsolve flag combos
nnz_combos = cat(1, ident_solution.flag)~=0;
nnz_combo_id = find(nnz_combos);
% collect all nnz parameter values
p_val = cat(2, ident_solution(nnz_combos).p);
[np, n_nnz_combos] = size(p_val);
p_val_signs = sign(p_val);
pos_combos = zeros(np,n_nnz_combos);
for i_p = 1:np
    pos_combos(i_p, p_val_signs(i_p, :)>0) = 1;    
end

nnz_flags = cat(1, ident_solution(nnz_combos).flag);
% get data for fully solved combos
nnz_sol = nnz_flags==1;
nnz_sol_id = nnz_combo_id(nnz_sol);
p_val_nnz_sol = p_val(:, nnz_sol);
pos_combos_nnz_sol = pos_combos(:, nnz_sol);

data = struct();
data.full_solved.combo_id = nnz_sol_id;
data.full_solved.p = p_val_nnz_sol;
data.full_solved.pos_p_combos = pos_combos_nnz_sol;


% get box plots
% figure
% boxplot(p_val')

% collect all solved combinations: flag == 2 TolFun < specified TolFun
% solved_combos = cat(1, ident_solution.flag)==2;
% collect all solved combinations: flag == 3 inaccurancy possible
% inaccurate_solve_combos = cat(1, ident_solution.flag==3);
% collect all unsolved combinations: flag == 0
% unsolved_combos = cat(1, ident_solution.flag)==0;
% collect all unsolved combinations: flag == -2
% ineffect_last_step_combos = cat(1, ident_solution.flag)==-2;
% collect all unsolved combinations: flag == -3
% singular_unsolved_combos = cat(1, ident_solution.flag)==-3;

% get results from fully solved combos
% full_solved_values = cat(2, ident_solution(full_solved_combos).p);


return