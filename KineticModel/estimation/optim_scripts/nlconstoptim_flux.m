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
if ~isempty(objCASh)
    c = sparse(7,1,-1,7,1); % [0 0 0 0 0 0 -1];
    [FXobj,gradF] = objCASh();
    obj = @(x)full(FXobj(x,c));
    grad = @(x)full(gradF(x,c));
end

% constraint function
if ~isempty(constrfh)
    % nlcons = @(x)constrfh(x,optim_p);
    [FXcons,DFXcons] = constrfh();
    nlconFX = @(x)full(FXcons(x,optim_p));
    nlconDFX = @(x)full(DFXcons(x,optim_p));
    givenconstr = 1;
else
    givenconstr = 0;
end

% constraint rhs
nlrhs = zeros(m+2,1);
nlrhs(1:2) = 0.1;
% nlrhs(3:6) = 1e-6; % precision level for Sv <= 0
% constraint type : (-1 <=, 0 =, 1 >=)
nle = zeros(m+2,1); % -ones(m+2,1);
nle(1:2) = -1;

% test if x0 satisfies constraints 
nlconsval = nlconFX(x0);
% == constraints
eqcons_iid = ones(length(nlrhs),1);
if ~all(nlconsval(nle==0) == nlrhs(nle==0))
    eqcons_iid(nlconsval ~= nlrhs) = -1;
    eqcons_iid(nle~=0) = 0;    
end
% <= constraints
lecons_iid = ones(length(nlrhs),1);
if ~all(nlconsval(nle==-1) == nlrhs(nle==-1))
    lecons_iid(nlconsval > nlrhs) = -1;
    lecons_iid(nle~=-1) = 0;    
end
% >= constraints
gecons_iid = ones(length(nlrhs),1);
if ~all(nlconsval(nle==1) == nlrhs(nle==1))
    gecons_iid(nlconsval < nlrhs) = -1;
    gecons_iid(nle~=1) = 0;    
end

% if any(eqcons_iid(nle==0)<0) ||...
%     any(lecons_iid(nle==-1)<0) ||...
%     any(gecons_iid(nle==1)<0)    
%     error('Initial point is infeasible - constraint(s) are violated');
% end
        


% setup opti problem
solver_opts = nloptset('algorithm','AUGLAG','subalgorithm','LN_COBYLA'); 
optim_opts = optiset('solver','NLOPT','maxiter',5000,...
                                      'maxfeval',500000,...
                                      'tolrfun',1e-8,...
                                      'tolafun',1e-10,...
                                      'solverOpts',solver_opts,...
                                      'display','final');
% solver_opts
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


