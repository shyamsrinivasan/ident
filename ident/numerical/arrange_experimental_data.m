% arrange experimental data (parameters) to perform numerical analysis of
% identifiability for v3 in small kotte model
function ident_data = arrange_experimental_data(perturbation_data,...
                                                nexpts,...
                                                select_exp_id)
if nargin<3
    select_exp_id = [];
end
nperts = size(perturbation_data, 2);
% generate combinations of experiments using combinator from matlab central
% combinations are not ordered
all_combos = combinator(nperts, nexpts, 'p');

% select only combos w/ or w/o specific experiment ids
if ~isempty(select_exp_id)
    chosen_combos = select_combos(all_combos, nexpts, select_exp_id);
end
return


function combos_w_exp_id = select_combos(all_combos, nexpts, select_exp_id)
% select only combos that have select_exp_id in them
useful_combos = zeros(size(all_combos, 1), length(select_exp_id));
for i_exp_id = 1:length(select_exp_id)
    combos_w_exp_id = [];
    for i_position = 1:nexpts
        combos_w_exp_id = union(combos_w_exp_id,...
                                find(all_combos(:, i_position)));
    end
    useful_combos(combos_w_exp_id, i_exp_id) = 1;    
end
return


