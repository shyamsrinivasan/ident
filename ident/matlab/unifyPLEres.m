% uify results from prallel computation of PLE
function newPLE = unifyPLEres(oldPLE)
% oldPLE has 2 cells - +ve and -ve perturbations

if ~isempty(oldPLE{1}) && ~isempty(oldPLE{2})
    field_names = fieldnames(oldPLE{1});
    for jnames = 1:length(field_names)
        f_name = field_names{jnames};
        if ~strcmpi(f_name,'iter')&&~strcmpi(f_name,'pos_info')&&~strcmpi(f_name,'neg_info')
            newPLE.(f_name) = [oldPLE{2}.(f_name) oldPLE{1}.(f_name)];
        elseif strcmpi(f_name,'iter')
            newPLE.(f_name) = [oldPLE{1}.(f_name)(1) oldPLE{2}.(f_name)(2)];
        elseif strcmpi(f_name,'pos_info')
            if ~isempty(oldPLE{1}.(f_name))
                newPLE.(f_name) = oldPLE{1}.(f_name);
            else
                newPLE.(f_name) = [];
            end            
        elseif strcmpi(f_name,'neg_info')
            if ~isempty(oldPLE{2}.(f_name))
                newPLE.(f_name) = oldPLE{2}.(f_name);
            else
                newPLE.(f_name) = [];
            end
        end
    end        
else
    newPLE =...
    oldPLE{~logical(cell2mat((cellfun(@isempty,oldPLE,'UniformOutput',false))))};
end



