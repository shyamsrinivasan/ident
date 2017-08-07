% estimate fluxes for flux 2 in kotte model pep ---> fdp

% generate experimental data - get initial ss
tspan = 0:0.1:100;
[xdyn,fdyn,xss,fss,opts] = run_nonoise(tspan);

% backup parameters and initial conditions
ival_bkp = opts.x0;
odep_bkp = opts.odep;

% perturb system from non-noisy initial conditions
pt_sol_id = [1 2 3];
[exp_sol,no_noise_sol] = dopert_nonoise(opts,xss,fss,odep_bkp,pt_sol_id);
close all

% get only data from one steady state
pss = ones(1,numel(exp_sol.exp_pval));
% pss(exp_sol.xss(1,:)>exp_sol.xss(2,:)) = 0;    

% options structure for solving problem
optimdata = struct('nc',3,'nflx',6,'nf',1,'flxid',4,'eps_v',1,...
                    'eps_c',1,'vexp',exp_sol.fss(:,logical(pss)),...
                    'xexp',exp_sol.xss(:,logical(pss)),...
                    'flux_wt',1000,'conc_wt',100,...
                    'eps_c_wt',1000,'eps_v_wt',10000,...
                    'odep',odep_bkp,...
                    'wt_xss',xss(:,1),'wt_fss',fss(:,1),...
                    'type',2);
                
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
% define all constraint fun handles consnoisyf2, cons_typeb_f2
allconsh = {[],...
            [],...
            [],...
            'cons_typeb_f2',...
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

[prob,optimdata] = setup_optim_prob(setup_opts);

% initial values for consrained nl(or quadratic?) optimization
x0 = getrandomivals(optimdata,.3,1000);
solveropt = struct('solver','ipopt','multi',0);
optsol = choose_nlconsopt(prob,x0,optimdata,solveropt);

% combine results for comparison plot
opts.tspan = 1:.1:200;  
est_data = combine_results(optsol,opts,no_noise_sol,optimdata,pss,pss);

% compare fluxes and concentrations
hfcv = compare_vals(est_data,no_noise_sol,optimdata,opts,pss);

% compare parameters in parameter space
hfp = compare_pars(est_data);

% save figure files
dir = 'C:\Users\shyam\Documents\Courses\CHE1125Project\Results\estimation\est_flux2\no_noise\typec\';
if ~isempty(hfcv)
    set(0,'CurrentFigure',hfcv(1));
    fname = 'est_flux2_conc_aug6';
    print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
    close(hfcv(1));
    set(0,'CurrentFigure',hfcv(2));
    fname = 'est_flux2_flux_aug6';
    print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
    close(hfcv(2));
end
if ~isempty(hfp)
    set(0,'CurrentFigure',hfp);
    fname = 'est_flux2_par_aug6';
    print([dir fname],'-depsc','-painters','-loose','-tiff','-r200');
    close(hfp);
end




