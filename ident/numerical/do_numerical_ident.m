function solution = do_numerical_ident(exp_data,ident_struct)
ndata = size(exp_data,1);
% get nlae solution for every data set in nadata
% consider paralleling
if isfield(ident_struct,'fun')
    ident_fun = ident_struct.fun;
end
if isfield(ident_struct,'initial_value')
    initial_val = ident_struct.initial_value;
end
if isfield(ident_struct,'typical_value')
    typical_val = ident_struct.typical_value;
end

solution = struct([]);
% assign structure fields in case of parallelization
for idata = 1:ndata
    solution(idata).p = [];
    solution(idata).res = [];
    solution(idata).flag = [];
    solution(idata).message = [];
end
parfor idata = 1:ndata
    single_solution =...
    nlae_solution(ident_fun,exp_data(idata, :),initial_val, typical_val);
    solution(idata).p = single_solution.p;
    solution(idata).res = single_solution.res;
    solution(idata).flag = single_solution.flag;
    solution(idata).message = single_solution.message;
    fprintf('Completed identification for data %d of %d', idata, ndata);
end

return