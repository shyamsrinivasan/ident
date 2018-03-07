function hf = get_boxplots(data)
% get info on positive parameter values only for each parameter
all_p_value = data.p;
pos_p_value_id = data.pos_p_combos;
np = size(all_p_value, 1);
pos_p_value = cell(np, 1);
pos_p_group = cell(np, 1);
for ip = 1:np
    parameter_name = sprintf('parameter %d',ip);
    pos_p_value{ip, 1} = all_p_value(ip, logical(pos_p_value_id(ip, :)));
    pos_p_group{ip, 1} =...
        repmat({parameter_name}, 1, length(pos_p_value{ip, 1}));
end
all_data = [pos_p_value{:}]';
all_group = [pos_p_group{:}]';

hf = figure;
boxplot(all_data, all_group);
return