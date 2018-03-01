% arrange experimental data (parameters) to perform numerical analysis of
% identifiability for v3 in small kotte model
function chosen_data = arrange_experimental_data(perturbation_data,...
                                                nexpts,...
                                                flux_id,...
                                                select_exp_id)
if nargin<4
    select_exp_id = [];
end
nperts = size(perturbation_data, 2);
% generate combinations of experiments using combinator from matlab central
% combinations are not ordered
all_combos = combinator(nperts, nexpts, 'p');

% select only combos w/ or w/o specific experiment ids
if ~isempty(select_exp_id)
    chosen_combos = select_combos(all_combos, nexpts, select_exp_id);
else
    chosen_combos = all_combos;
end

% get data for chosen combos
chosen_data = choose_experimental_data(chosen_combos,perturbation_data,nexpts,flux_id);
return

function selected_combos = select_combos(all_combos, nexpts, select_exp_id)
% select only combos that have select_exp_id in them
useful_combos = zeros(size(all_combos, 1), length(select_exp_id));
for i_exp_id = 1:length(select_exp_id)
    combos_w_exp_id = [];
    for i_position = 1:nexpts
        combos_w_exp_id =...
            union(combos_w_exp_id,...
                  find(all_combos(:, i_position)==select_exp_id(i_exp_id)));
    end
    useful_combos(combos_w_exp_id, i_exp_id) = 1;     
end
useful_combo_id = logical(sum(useful_combos, 2));
selected_combos = all_combos(useful_combo_id, :);
return

function chosen_data = choose_experimental_data(chosen_combos,...
                                                perturbation_data,...
                                                nexpts,...
                                                flux_id)
% select data from experiments in chosen_combos for each combination in 
% chosen_combos
% 3 concentrations(pep, fdp, E) + 4 fluxes(v1,v2,v3,v4)+ 1
% parameter(acetate)
data_size = 8; 
ncombos = size(chosen_combos,1);
chosen_data = zeros(ncombos, data_size*nexpts);
for icombo = 1:ncombos
    chosen_data(icombo, :) =...
        select_experimental_data(perturbation_data,...
                                 chosen_combos(icombo, :),...
                                 flux_id);    
end
return

function dataset =...
         select_experimental_data(perturbation_data, exp_id, flux_id)
% get data for each experiment within each chosen combination passed as
% input
if nargin<3
    flux_id = 0:length(perturbation_data(1).fss);
end
dataset = struct();
nexpts = length(exp_id);
for iexp = 1:nexpts
    dataset(iexp).data = [perturbation_data(exp_id(iexp)).parameter(end);...
                          perturbation_data(exp_id(iexp)).yss;...
                          perturbation_data(exp_id(iexp)).fss(flux_id)];
end
dataset = cat(1, dataset.data);
return



