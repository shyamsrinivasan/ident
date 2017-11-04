% use perturbations in xss and fss to get flux parameters
function [optsol] = optimize_p_noisy(opts,xss,fss,plist,odep_bkp)
% if isfield(opts,'init_xss')
%     init_xss = opts.init_xss;
% end
opts.odep = odep_bkp;
optsol = struct([]); % cell(4,3);
% opt_id = cell(4,1);

ts1 = tic;

% 'K1ac','k1cat'
fprintf('Optimizing with a single set of perturbations\n');
fprintf('\nOptimizing parameters for flux 1.....\n');
[x_opt_4_1,opt_id_4_1,~,fval_4_1] = flux1_k_noisy(opts,xss,fss,plist);
optsol(1).x_opt = x_opt_4_1;
optsol(1).opt_id = opt_id_4_1;
optsol(1).fval = fval_4_1;
optsol(1).time = toc(ts1);
fprintf('Time to estimate flux 1 parameters %4.3g\n',toc(ts1));

% % 'K2pep','V2max','K2fdp'
% ts2 = tic;
% fprintf('\nOptimizing parameters for flux 2.....\n');
% [x_opt_4_2,opt_id_4_2,~,fval_4_2] = flux2_k(opts,xss,fss,plist);
% optsol(2).x_opt = x_opt_4_2;
% optsol(2).opt_id = opt_id_4_2;
% optsol(2).fval = fval_4_2;
% optsol(2).time = toc(ts2);
% fprintf('Time to estimate flux 2 parameters %4.3g\n',toc(ts2));
% 
% % 'K3fdp','K3pep','V3max'
% ts3 = tic;
% fprintf('\nOptimizing parameters for flux 3.....\n'); 
% [x_opt_4_3,opt_id_4_3,~,fval_4_3] = flux3_k(opts,xss,fss,plist);
% optsol(3).x_opt = x_opt_4_3;
% optsol(3).opt_id = opt_id_4_3;
% optsol(3).fval = fval_4_3;
% optsol(3).time = toc(ts3);
% fprintf('Time to estimate flux 3 parameters %4.3g\n',toc(ts3));
% 
% % 'V4max'
% ts4 = tic;
% fprintf('\nOptimizing parameters for flux 4.....\n');
% [x_opt_4_4,opt_id_4_4,~,fval_4_4] = flux4_k(opts,xss,fss,plist);
% optsol(4).x_opt = x_opt_4_4;
% optsol(4).opt_id = opt_id_4_4;
% optsol(4).fval = fval_4_4;
% optsol(4).time = toc(ts4);
% fprintf('Time to estimate flux 1 parameters %4.3g\n',toc(ts4));

fprintf('\nTotal time for parameter estimation %4.3g\n',toc(ts1));