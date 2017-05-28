function [x_opt,fval,xss4,fss4,opts] =...
         runoptim_flux(opts,objCASh,lb,ub,x0,optim_p,multi)

if nargin<7
    multi = 0;
end

[FX,gradF] = objCASh();
obj = @(x)full(FX(x,optim_p));
grad = @(x)full(gradF(x,optim_p));

% setup opti problem
solver_opts = nloptset('algorithm','GD_STOGO_RAND');
optim_opts = optiset('solver','NLOPT','maxiter',5000,...
                                      'maxfeval',5000,...
                                      'tolrfun',1e-8,...
                                      'tolafun',1e-10,...
                                      'solverOpts',solver_opts,...
                                      'display','final');
% optim_opts = optiset('solver','IPOPT','maxiter',5000,...
%                                       'maxfeval',5000,...                                      
%                                       'display','final');
optim_prob = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',optim_opts);
if multi
    [xval,fval,exitflag,info] = multisolve(optim_prob,[],[50 10]);   
else
    [xval,fval,exitflag,info] = solve(optim_prob,x0); 
end

if exitflag == 1 &&...
   (strcmpi(info.Status,'Success')||...
    strcmpi(info.Status,'Optimal')||...
    strcmpi(info.Status,'Converged / Target Reached'))

    x_opt = xval;
else
    x_opt = [];
    fval = []; 
end

xss4 = [];
fss4 = [];

