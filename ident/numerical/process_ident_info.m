% process numerical ientifiability information after solving systems of
% nlae to identify parameters
function process_ident_info(ident_solution)
% collect all solved combinations
full_solved_combos = cat(1, ident_solution.flag)==1;
% collect all solved combinations: flag == 2 TolFun < specified TolFun
solved_combos = cat(1, ident_solution.flag)==2;
% collect all solved combinations: flag == 3 inaccurancy possible
inaccurate_solve_combos = cat(1, ident_solution.flag==3);
% collect all unsolved combinations: flag == 0
unsolved_combos = cat(1, ident_solution.flag)==0;
% collect all unsolved combinations: flag == -2
ineffect_last_step_combos = cat(1, ident_solution.flag)==-2;
% collect all unsolved combinations: flag == -3
singular_unsolved_combos = cat(1, ident_solution.flag)==-3;

% get results from fully solved combos
full_solved_values = cat(2, ident_solution(full_solved_combos).p);


return