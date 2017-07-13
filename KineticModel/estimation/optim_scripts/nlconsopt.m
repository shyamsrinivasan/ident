function optsol = nlconsopt(prob,x0,solveropts,miscdata)
if nargin<4
    miscdata = [];
end
% problem structure
if isfield(prob,'obj')
    obj = prob.obj;
else
    obj = [];
end
if isfield(prob,'nlcons')
    nlcons = prob.nlcons;
else
    nlcons = [];
end
if isfield(prob,'lb')
    lb = prob.lb;
else
    lb = [];
end
if isfield(prob,'ub')
    ub = prob.ub;
else
    ub = [];
end
% if isfield(prob,'optimopts')
%     optimopts = prob.optimopts;
% else
%     optimopts = [];
% end
if isfield(prob,'nlrhs')
    nlrhs = prob.nlrhs;
else
    nlrhs = [];
end
if isfield(prob,'nle')
    nle = prob.nle;
else
    nle = [];
end

% options structure
if isfield(solveropts,'solver')
    solver = solveropts.solver;
else
    solver = 'nlopt';
end
if isfield(solveropts,'opt_alg')
    opt_alg = solveropts.opt_alg;
else
    opt_alg = 'none';
end
if isfield(solveropts,'multi')
    multi = solveropts.multi;
else
    multi = 0;
end

if isempty(obj)
    error('No given objective function');
end

tstrt = tic;
optimopts = optiset('solver',solver,...
                      'maxiter',5000,...
                      'maxfeval',500000,...
                      'tolrfun',1e-6,...
                      'tolafun',1e-6,...
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
    case 'ipopt'
         solver_opts = [];
    case 'scip'
        solver_opts = scipset('scipopts',{'limits/time',1e6});
        optimopts = optiset(optimopts,'maxnodes',100000000);
    otherwise
        error('Nonexistent solver in OPTI');
end
optimopts = optiset(optimopts,'solverOpts',solver_opts);   

if ~isempty(nlcons)
    optim_prob =...
    opti('obj',obj,'nlmix',nlcons,nlrhs,nle,...
         'ndec',miscdata.nvar,'bounds',lb,ub,'options',optimopts);
else
    optim_prob =...
    opti('obj',obj,'bounds',lb,ub,'options',optimopts);    
end
if multi
    [xval,fval,exitflag,info] = multisolve(optim_prob,[],[3 10]);   
else
    [xval,fval,exitflag,info] = solve(optim_prob,x0); 
end     
% [xval,fval,exitflag,info] = solve(optim_prob,x0); 

optsol.x0 = x0;
optsol.xval = xval;
optsol.fval = fval;
optsol.exitflag = exitflag;
optsol.info = info;
optsol.time = toc(tstrt);