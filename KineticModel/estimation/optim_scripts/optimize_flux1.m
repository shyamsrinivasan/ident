function optsol = optimize_flux1(opts,xss,fss,plist,odep_bkp,pid,pval)
if nargin<7
    pval = [];
end
if nargin<6
    pid = [];
end

opts.odep = odep_bkp;
optsol = struct([]); 

ts1 = tic;

% 'K1ac','k1cat'
fprintf('Optimizing with a single set of perturbations\n');
fprintf('\nOptimizing parameters for flux 1.....\n');

% flux 1    
p_id = cellfun(@(x)strcmpi(plist,x),{'K1ac','k1cat'},'UniformOutput',false);
p_id = cellfun(@(x)find(x),p_id);

% set initial values x0
% flux 1 estimation
% p_id = [1,11];               
% ival_id = 1:4;
nss_var = 1:4;
% get average values of all experimental observations used for parameter
% estimation

xi = sum(xss(:,nss_var),2)/length(nss_var);
pi = opts.odep;
pi(pid) = sum(pval(nss_var))/length(nss_var);
x0 = [xi;pi(p_id)';0];
opts.opt_x0 = x0;

ss_val = [xss(:,nss_var);...
          fss(1,nss_var)];
      
% check feasibility of x0
% cons_rhs = constr_flux1_noisy(x0,opts.odep,p_id,ss_val);


% optimization
[x_opt_4_1,opt_id_4_1,~,fval_4_1] = flux1_k_noisy(opts,ss_val,p_id);
optsol(1).x_opt = x_opt_4_1;
optsol(1).opt_id = opt_id_4_1;
optsol(1).fval = fval_4_1;
optsol(1).time = toc(ts1);
optsol(1).x0 = opts.opt_x0;
fprintf('Time to estimate flux 1 parameters %4.3g\n',toc(ts1));

fprintf('\nTotal time for parameter estimation %4.3g\n',toc(ts1));
