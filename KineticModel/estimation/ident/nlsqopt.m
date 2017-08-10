% optimization step for PLE-based identifiability
% idx - index of parameter that is fixed
function optsol = nlsqopt(prob,solveropts)
if isfield(prob,'x0')
    x0 = prob.x0;
end
if isfield(prob,'obj')
    obj = prob.obj;
else
    obj = [];
end
if isfield(prob,'gradobj')
    gradobj = prob.gradobj;
else
    gradobj = [];
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
if isfield(solveropts,'multi_pts')
    multi_pts = solveropts.multi_pts;
else
    multi_pts = [5 12];
end
if isempty(obj)
    error('No given objective function');
end

tstrt = tic;
optimopts = optiset('solver',solver,...
                      'maxiter',10000,...
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
    opti('fun',obj,'f',gradobj,'nlmix',nlcons,nlrhs,nle,...
         'ndec',miscdata.nvar,'bounds',lb,ub,'options',optimopts);
else
    optim_prob =...
    opti('fun',obj,'f',gradobj,'bounds',lb,ub,'options',optimopts);    
end
if multi
    [xval,fval,exitflag,info] = multisolve(optim_prob,[],multi_pts);   
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

% xopt = optsol.xval;

% info
% if isfield(data,'xexp')
%     yexp = data.xexp;
% end
% if isifield(data,'xexp_var')
%     yexp_var = data.yexp_var;
% end
% create nlsqopt objective from modelh (modelh is a casadi fun 
% that gives a time course profile of y = f(y,p))
% ymodel = modelh(x,p);
% obj = sum((yexp-ymodel)./yexp_var).^2;


% using casadi to solve nonlinear (objective) optimization 
% objective : y* - y
% bounds : p
% nlp = struct('x',var,'f',obj,'g',cons);
% solver = nlpsol('sol','ipopt',nlp);