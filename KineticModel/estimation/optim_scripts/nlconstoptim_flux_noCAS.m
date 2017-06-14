function [x_opt,fval,xss4,fss4,opts] =...
         nlconstoptim_flux_noCAS(opts,objCASh,lb,ub,x0,optim_p,ss_val,multi,constrfh)
if nargin<9
    givenconstr = 0;
else
    givenconstr = 1;
end
if nargin<8
    multi = 0;
end
if isfield(opts,'solver')
    solver = opts.solver;
else
    solver = 'nlopt';
end
if isfield(opts,'opt_alg')
    opt_alg = opts.opt_alg;
else
    opt_alg = 'none';
end
if isfield(opts,'nlrhs')
    nlrhs = opts.nlrhs;
end
if isfield(opts,'nle')
    nle = opts.nle;
end
np = size(ss_val,2); % # perturbations
m = 3;
n = 6;
nvar = 6;
var_id = 6;

% objective fun
if ~isempty(objCASh)
    c = sparse(var_id,1,1,nvar,1); % [0 0 0 0 0 0 1];
    obj = @(x)c'*x;
%     obj = @(x)objCASh(x,c);
%     grad = @(x)full(gradF(x,c));
end

% constraint function
if givenconstr
    nlconFX = @(x)constrfh(x,optim_p,opts.p_id,ss_val);     
end

% test if x0 satisfies constraints 
[eqcons_iid,lecons_iid,gecons_idd] = checkconstr(nlconFX,x0,nlrhs,nle);
% if any(eqcons_iid(nle==0)<0) ||...
%     any(lecons_iid(nle==-1)<0) ||...
%     any(gecons_iid(nle==1)<0)    
%     error('Initial point is infeasible - constraint(s) are violated');
% end

% setup opti problem
optim_opts = optiset('solver',solver,...
                      'maxiter',5000,...
                      'maxfeval',500000,...
                      'tolrfun',1e-8,...
                      'tolafun',1e-8,...
                      'display','final');                               
switch(solver)
    case 'nlopt'        
        switch(opt_alg)            
            case 'GN_ISRES'
                solver_opts = nloptset('algorithm','GN_ISRES'); 
            otherwise
                solver_opts =...
                nloptset('algorithm','AUGLAG','subalgorithm','LN_COBYLA'); 
        end
        optim_opts = optiset(optim_opts,'solverOpts',solver_opts);           
    case 'ipopt'
         
    case 'scip'
        optim_opts = optiset(optim_opts,'maxnodes',10000000);
    otherwise
        error('Nonexistent solver in OPTI');
end

if givenconstr
    optim_prob =...
    opti('obj',obj,'nlmix',nlconFX,nlrhs,nle,...
         'ndec',nvar,'bounds',lb,ub,'options',optim_opts);
else
    optim_prob =...
    opti('obj',obj,'bounds',lb,ub,'options',optim_opts);
end
if multi
    [xval,fval,exitflag,info] = multisolve(optim_prob,[],[5 10]);   
else
    [xval,fval,exitflag,info] = solve(optim_prob,x0); 
end     


if exitflag == 1 &&...
   (strcmpi(info.Status,'Success')||...
    strcmpi(info.Status,'Optimal')||...
    strcmpi(info.Status,'Converged / Target Reached')||...
    strcmpi(info.Status,'Globally Optimal'))

    x_opt = xval;
elseif ~isempty(xval) && ~isnan(fval)
    x_opt = xval;    
else
    x_opt = [];    
end

xss4 = [];
fss4 = [];


