function [x_opt,fval,xss4,fss4,opts] =...
        sciptest(opts,objh,lb,ub,x0,optim_p)
    
obj = @(x)objh(x,optim_p);    
% setup opti problem
solver_opts = nloptset('algorithm','GN_DIRECT'); % 'solver','pswarm',
optim_opts = optiset('solver','nlopt','maxiter',200000000,...
                                      'maxfeval',500000000,...
                                      'tolrfun',1e-10,...
                                      'tolafun',1e-12,...    
                                      'solverOpts',solver_opts,...
                                      'display','iter');
optim_prob = opti('obj',obj,'bounds',lb,ub,'options',optim_opts);
% [xval,fval,exitflag,info] = solve(optim_prob,x0);    
[xval,fval,exitflag,info] = multisolve(optim_prob,[],[50 10]);    

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


