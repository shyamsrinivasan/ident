% combine results from multiple optimization problems run for multi start
% search for comparison
function combine_results(optsol,exp_data,data,given_data_id,test_data_id)
% taken from parts of compare_vals
% given_data_id - sets in exp_data used for estimation
% test_data_id - sets in exp_data used for testing estimated parameters
if nargin<5
    test_data_id = given_data_id;
end

if ~isfield(optsol,'xconc')
    optsol = parse_optsol(optsol,data);
end

nval = size(optsol,2);
for ival = 1:nval
    est_xss
end


 