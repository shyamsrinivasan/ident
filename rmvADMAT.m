function rmvADMAT(tool)

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
if any(ind)
    fprintf('path %s is in MATLAB search path\n',tool_path);
    fprintf('All folders and subfolders of \n %s will be removed\n',tool_path);
    fprintf('Removing paths...\n');
    disp(path_list(ind));
    rmpath(path_list(ind));
end