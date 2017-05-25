function [x_opt,fval,xss4,fss4,opts] =...
         runoptim_flux(opts,objCASh,lb,ub,x0,optim_p)

[FX,gradF] = objCASh();
obj = @(x)full(FX(x,optim_p));
grad = @(x)full(gradF(x,optim_p));

% setup opti problem
optim_opts = optiset('solver','NLOPT','maxiter',5000,'maxfeval',5000,'display','iter');
optim_prob = opti('obj',obj,'grad',grad,'bounds',lb,ub,'options',optim_opts);
[xval,fval,exitflag,info] = solve(optim_prob,x0); 

if exitflag == 1 &&...
   (strcmpi(info.Status,'Success')||strcmpi(info.Status,'Optimal'))
    x_opt = xval;
else
    x_opt = [];
    fval = []; 
end

xss4 = [];
fss4 = [];

