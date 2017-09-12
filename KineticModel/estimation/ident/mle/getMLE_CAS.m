% get maximum likelihood estimate(mle) of parameters to begin profile
% likelihood analysis
% log L = -1/2SUM(yi - yhati)^2/sigma^2
% MLE = min log(L)
function MLEvals = getMLE_CAS(data,x0,scale)

prob_struct = mle_setupCAS(data);
optsol = solve_nlsqopt(prob_struct,x0);

optsol.xval = optsol.xval.*scale;
mle_pval = [optsol.xval(1:5);data.odep(6:9)';optsol.xval(6:9);data.odep(14:17)'];
MLEvals.mle_pval = mle_pval;
MLEvals.logL = optsol.fval;
% MLEvals.exitflag = exitflag;
MLEvals.info = optsol.info;



