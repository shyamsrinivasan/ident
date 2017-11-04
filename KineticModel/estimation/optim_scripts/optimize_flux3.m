function optsol = optimize_flux3(opts,xss,fss,plist,odep_bkp)

opts.odep = odep_bkp;
optsol = struct([]); 

ts1 = tic;

% 'K1ac','k1cat'
fprintf('Optimizing with a single set of perturbations\n');
fprintf('\nOptimizing parameters for flux 3.....\n');

% set initial values x0
% flux 1 estimation
p_id = [2,4,12];               
xi = xss(:,end);
x0 = [xi;opts.odep(p_id)';.1];
opts.opt_x0 = x0;

% optimization
[x_opt_4_1,opt_id_4_1,~,fval_4_1] = flux3_k_noisy(opts,xss,fss,plist);
optsol(1).x_opt = x_opt_4_1;
optsol(1).opt_id = opt_id_4_1;
optsol(1).fval = fval_4_1;
optsol(1).time = toc(ts1);
fprintf('Time to estimate flux 2 parameters %4.3g\n',toc(ts1));

fprintf('\nTotal time for parameter estimation %4.3g\n',toc(ts1));
