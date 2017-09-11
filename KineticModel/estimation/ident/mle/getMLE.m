% get maximum likelihood estimate(mle) of parameters to begin profile
% likelihood analysis
% log L = -1/2SUM(yi - yhati)^2/sigma^2
% MLE = min log(L)
function getMLE(data,x0)

prob_struct = mle_setup(data);

% run optimization with x0 as initial parameter values

