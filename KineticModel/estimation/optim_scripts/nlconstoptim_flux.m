function [x_opt,fval,xss4,fss4,opts] =...
         nlconstoptim_flux(opts,objCASh,lb,ub,x0,optim_p,multi,constrfh)

if nargin<7
    multi = 0;
end
if nargin<8
    constrfh = [];
end
% np = size(optim_p,2);
m = 4;
% objective fun
obj = @(x)sparse(7,1,-1,7,1)'*x;
grad = @(x)-1;
% c = ; % [-1 0 0 0 0];

% constraint function
if ~isempty(constrfh)
    % nlcons = @(x)constrfh(x,optim_p);
    [FX,DFX] = constrfh();
    nlconFX = @(x)full(FX(x,optim_p));
    nlconDFX = @(x)full(DFX(x,optim_p));
    givenconstr = 1;
else
    givenconstr = 0;
end
% test if x0 satisfies constraints constraint fun
nlconsval = nlconFX(x0);
% constraint rhs
nlrhs = zeros(m+2,1);
nlrhs(1:2) = 0.1;
nlrhs(3:6) = 1e-10; % precision level for Sv <= 0
% constraint type : (-1 <=, 0 =, 1 >=)
nle = -ones(m+2,1); % zeros(m+2,1);
% nle(1:2) = -1;

% FX = objCASh(np);
% obj = @(x)full(FX(x,optim_p));
% grad = @(x)full(gradF(x,optim_p));

% setup opti problem
% solver_opts = nloptset('algorithm','GN_ISRES');'solverOpts',solver_opts,...
optim_opts = optiset('solver','NLOPT','maxiter',5000,...
                                      'maxfeval',500000,...
                                      'tolrfun',1e-8,...
                                      'tolafun',1e-10,...                                      
                                      'display','final');
% optim_opts = optiset('solver','IPOPT','maxiter',5000,...
%                                       'maxfeval',5000,...                                      
%                                       'display','final');
if givenconstr
    optim_prob =...
    opti('obj',obj,'nlmix',nlconFX,nlrhs,nle,'nljac',nlconDFX,...
         'ndec',7,'bounds',lb,ub,'options',optim_opts);
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


