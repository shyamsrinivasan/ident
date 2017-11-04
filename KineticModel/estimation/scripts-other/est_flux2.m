% estimate fluxes for flux 2 in kotte model pep ---> fdp

% run gen_data for experimental data
%% gen_data

% or load pre-calculated data
load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noisy_model/pdata_aug24');

%% run optimization
id = 1;
% get only data from one steady state
pss = ones(1,numel(exp_sol{id}.exp_pval));
% pss(exp_sol.xss(1,:)>exp_sol.xss(2,:)) = 0;    

% options structure for solving problem
optimdata = struct('nc',3,'nflx',6,'nf',1,'flxid',4,'eps_v',1,...
                    'eps_c',1,'vexp',exp_sol{id}.fss(:,logical(pss)),...
                    'xexp',exp_sol{id}.xss(:,logical(pss)),...
                    'flux_wt',1000,'conc_wt',1,...
                    'eps_c_wt',1,'eps_v_wt',1000,...
                    'odep',odep_bkp,...
                    'wt_xss',noisy_xss(:,1),'wt_fss',noisy_fss(:,1),...
                    'type',2);
expdata = struct('vexp',exp_sol{id}.fss(:,logical(pss)),...
                'xexp',exp_sol{id}.xss(:,logical(pss)));
                
% define all bound fun handles
allboundh = {[],...
             [],...
             [],...
             'boundsf2',...
             [],...
             []};
% define all objectiv efun handles
allobjh = {[],...
           [],...
           [],...
           'obj_typec',...
           [],...
           []};
% define all constraint fun handles
allconsh = {[],...
            [],...
            [],...
            'cons_typec_f2',...
            [],...
            []};
% define all rhs for nl cons
allnlrhs = {0,[],0,[],0,[]};
% define cons type for nl cons
allnle = {0,[],0,[],0,[]};

% problem defn
setup_opts = optimdata;
setup_opts.bounds = allboundh;
setup_opts.obj = allobjh;
setup_opts.nlcons = allconsh;
setup_opts.nlrhs = allnlrhs;
setup_opts.nle = allnle;

[prob,optimdata] = setup_optim_prob(setup_opts,expdata);

% initial values for consrained nl(or quadratic?) optimization
x0 = getrandomivals(optimdata,.3,5000);
solveropt = struct('solver','ipopt','multi',0);
optsol = choose_nlconsopt(prob,x0,optimdata,solveropt);

%% parse results
% load stored data
load('C:/Users/shyam/Documents/Courses/CHE1125Project/Results/estimation/mat_files/est_flux2_5000_aug22');

% combine results for comparison plot
opts.tspan = 1:.1:200;
[proc_data,noisy_sol] = recalcss(optsol,noisy_sol,[],optimdata,opts);
% est_data = combine_results(optsol,opts,noisy_sol,optimdata,pss,pss);

%% compare fluxes and concentrations
hfcv = compare_vals(proc_data,noisy_sol,[],optimdata,1);
hfdotcv = compare_vals(proc_data,noisy_sol,[],optimdata,2);

compare_vals_scatter(noisy_sol,[],proc_data.opt_xss,proc_data.calc_xss,optimdata,2);
compare_vals_scatter(noisy_sol,[],proc_data.opt_xss,proc_data.calc_xss,optimdata,1);

% compare parameters in parameter space
% hfp = compare_pars(est_data);

% save figure files
% dir = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\estimation\est_flux2\noisy\typec\';
% set(0,'CurrentFigure',hfcv(1));
% fname = 'est_flux2_conc_aug6';
% print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
% close(hfcv(1));
% set(0,'CurrentFigure',hfcv(2));
% fname = 'est_flux2_flux_aug6';
% print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
% close(hfcv(2));
% set(0,'CurrentFigure',hfp);
% fname = 'est_flux2_par_aug6';
% print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
% close(hfp);




