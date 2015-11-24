function addADMAT(tool)

base = 'C:\Users\shyam\Documents\MATLAB';
basepath = genpath(base);
basepath_list = regexp(basepath,';','Split');
basepath_list = basepath_list(~cellfun('isempty',basepath_list));

% tool = 'zz_ADMAT-2.0';
tool_path = strcat(base,'\',tool);

path = genpath(tool_path);
path_list = regexp(path,';','Split');
path_list = path_list(~cellfun('isempty',path_list));

ind = ismember(path_list,basepath_list);
if ~all(ind)
    fprintf('path %s not in MATLAB search path\n',tool_path);
    fprintf('All folders and subfolders of \n %s will be added\n',tool_path);
    fprintf('Adding paths...\n');
    disp(path_list(~ind));
    addpath(path_list(~ind),'-end');
end
    