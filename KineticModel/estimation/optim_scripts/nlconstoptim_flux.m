function [x_opt,fval,xss4,fss4,opts] =...
         nlconstoptim_flux(opts,objCASh,lb,ub,x0,optim_p,multi,constrfh)

if nargin<7
    multi = 0;
end
if nargin<8
    costrfh = [];
end
% np = size(optim_p,2);

% objective vector
e = sparse(1,1,-1,1,5); % [-1 0 0 0 0];

% constraint function
if ~isempty(costrfh)
    nlcons = @(x)constrfh(x,optim_p);
    givenconstr = 1;
else
    givenconstr = 0;
end
% constraint rhs
nlrhs
% constraint type : (-1 <, 0 =, 1 >)
nle

% FX = objCASh(np);
% obj = @(x)full(FX(x,optim_p));
% grad = @(x)full(gradF(x,optim_p));

% setup opti problem
solver_opts = nloptset('algorithm','GN_DIRECT');
optim_opts = optiset('solver','NLOPT','maxiter',5000,...
                                      'maxfeval',5000,...
                                      'tolrfun',1e-8,...
                                      'tolafun',1e-10,...
                                      'solverOpts',solver_opts,...
                                      'display','final');
% optim_opts = optiset('solver','IPOPT','maxiter',5000,...
%                                       'maxfeval',5000,...                                      
%                                       'display','final');
if givenconstr
    optim_prob = opti('obj',e,'nlmix',nlcons,nlrhs,nle,'bounds',lb,ub,'options',optim_opts);
else
    optim_prob =...
    opti('obj',obj,'bounds',lb,ub,'options',optim_opts);
end
if multi
    [xval,fval,exitflag,info] = multisolve(optim_prob,[],[100 10]);   
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


