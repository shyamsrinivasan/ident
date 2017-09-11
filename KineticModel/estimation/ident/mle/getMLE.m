% get maximum likelihood estimate(mle) of parameters to begin profile
% likelihood analysis
% log L = -1/2SUM(yi - yhati)^2/sigma^2
% MLE = min log(L)
function MLEvals = getMLE(data,x0,scale)

prob_struct = mle_setup(data);

% run optimization with x0 as initial parameter values
% solver_opts = ipoptset('max_cpu_time',1e8);
opti_options = optiset('solver','ipopt',...
              'maxiter',10000,...
              'maxfeval',500000,...
              'tolrfun',1e-6,...
              'tolafun',1e-6,...
              'display','iter',...
              'maxtime',3000);  
          
% use ipopt from opti instead of casadi version          
prob =...
opti('obj',prob_struct.objfun,...
     'grad',prob_struct.grad,...
     'bounds',prob_struct.lb,prob_struct.ub,...
     'options',opti_options);  
 
[xval,fval,exitflag,info] = solve(prob,x0);  

if exitflag==1
    p_opt = xval.*scale;
    logL = fval; % value of log likelihood
else
    p_opt = [];
    logL = fval;
end

if ~isempty(p_opt)
    % original parameter list in data.odep
%     plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep','vemax','KeFDP','ne',...
%             'd','V4max','k1cat','V3max','V2max','K1pep','K2fdp','rhoA','acetate'}; 
    % new parameter list for MLE
%     plist = {'K1ac','K3fdp','L3fdp','K3pep','K2pep',...
%             'V4max','k1cat','V3max','V2max'}; 
    mle_pval = [p_opt(1:5);data.odep(6:9)';p_opt(6:9);data.odep(14:17)'];
else
    mle_pval = [];
end
MLEvals.mle_pval = mle_pval;
MLEvals.logL = logL;
MLEvals.exitflag = exitflag;
MLEvals.info = info;