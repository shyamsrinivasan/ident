% combine results from multiple optimization problems run for multi start
% search for comparison
function est_data =...
combine_results(optsol,odeopts,exp_data,data,given_data_id,test_data_id)
% taken from parts of compare_vals
% given_data_id - sets in exp_data used for estimation
% test_data_id - sets in exp_data used for testing estimated parameters
if nargin<5
    test_data_id = given_data_id;
end
if ~islogical(given_data_id)
    given_data_id = logical(given_data_id);
end
if ~islogical(test_data_id)
    test_data_id = logical(test_data_id);
end

if ~isfield(optsol,'xconc')
    optsol = parse_optsol(optsol,data);
end

nval = size(optsol,2);
pt_datasize = length(find(test_data_id(1:end-1)));

exp_xss = cat(2,exp_data.xss);
exp_fss = cat(2,exp_data.fss);

wt_exp_xss = exp_xss(:,end);
wt_exp_fss = exp_fss(:,end);
% cell_optsol = cell(nval,1);

% assign same perturbations as required by test_data_id
[pt_val(1:pt_datasize).exp_pid] = exp_data(test_data_id(1:end-1)).exp_pid;
[pt_val(1:pt_datasize).exp_pval] = exp_data(test_data_id(1:end-1)).exp_pval;

% create nval matrix of parameters and nval structure array of options
odeopts_struct = struct();
opts_odep = repmat(data.odep,nval,1);
wt_est_xss_all = zeros(data.nc,nval);
for ival = 1:nval
    if ~test_data_id(end)
        wt_est_xss_all(:,ival) = wt_exp_xss; % model wt xss = exp wt xss         
    else
        wt_est_xss_all(:,ival) = optsol(ival).xconc(:,end);           
    end    
    
    % create options structure array
    opts_odep(ival,data.p_id) = optsol(ival).xpar;
    odeopts_struct(ival).tspan = odeopts.tspan;
    odeopts_struct(ival).x0 = odeopts.x0;
    odeopts_struct(ival).solver_opts = odeopts.solver_opts;
    odeopts_struct(ival).odep = opts_odep(ival,:);    
end

npert = data.npert;
parfor ival = 1:nval
    pss = ones(1,npert);
    % choose opts structure
    odeopts = odeopts_struct(ival);
    
%     if ~test_data_id(end)
%         wt_est_xss = wt_exp_xss; % model wt xss = exp wt xss         
%     else
%         wt_est_xss = optsol(ival).xconc(:,end);           
%     end    
%     xpar = optsol(ival).xpar;
%     opt_odep = data.odep;
%     opt_odep(data.p_id) = xpar;    
    
    wt_est_fss = kotte_flux_noCAS(wt_est_xss_all(:,ival),odeopts.odep);
    
%     odeopts.odep = opt_odep;
%     odeopts.tspan = 0:.1:500;
    odeopts.x0 = wt_est_xss_all(:,ival); 
    
    % recalculate perturbations for new parameters
    sol = getperturbations(pt_val,@perturb_nonoise,odeopts);
    close all
    
    % collect all solutions
%     cell_optsol{ival}.xss = [cat(2,sol.xss) wt_est_xss];
%     cell_optsol{ival}.fss = [cat(2,sol.fss) wt_est_fss];
    xss = [cat(2,sol.xss) wt_est_xss_all(:,ival)];
    fss = [cat(2,sol.fss) wt_est_fss];
    optsol(ival).xss = xss;
    optsol(ival).fss = fss;
    pss(xss(2,:)>xss(1,:)) = 2;
    optsol(ival).pss = pss;
end

% collect only optimal solutions
exitflags = cat(1,optsol.exitflag);
optimal_optsol = optsol(exitflags==1);

% recalculate nval based on optimal solutions only
nval = size(optimal_optsol,2);
npar = length(data.p_id);

% collect all concentrations and fluxes
conc = cat(1,optimal_optsol.xss);
flux = cat(1,optimal_optsol.fss);
par = cat(2,optimal_optsol.xpar);
est_xss = zeros(data.nc,data.npert);
est_xss_err = zeros(1,data.npert);
est_fss = zeros(data.nflx,data.npert);
est_fss_err = zeros(1,data.npert);
est_par = zeros(npar,1);
est_par_err = zeros(npar,1);
% concentrations
for ic = 1:data.nc
    est_xss(ic,:) = sum(conc(ic:data.nc:data.nc*nval,:),1)./nval;    
    est_xss_err(ic,:) = std(conc(ic:data.nc:data.nc*nval,:),0,1);
end
% fluxes
for iflx = 1:data.nflx
    est_fss(iflx,:) = sum(flux(iflx:data.nflx:data.nflx*nval,:),1)./nval;
    est_fss_err(iflx,:) = std(flux(iflx:data.nflx:data.nflx*nval,:),0,1);
end
% parameters
for ipar = 1:npar
    est_par(ipar) = sum(par(ipar,:))./nval;
    est_par_err(ipar) = std(par(ipar,:));
end
est_data.xconc = conc;
est_data.xss = est_xss;
est_data.xerr = est_xss_err;
est_data.xflux = flux;
est_data.fss = est_fss;
est_data.ferr = est_fss_err;
est_data.xpar = par;
est_data.par = est_par;
est_data.perr = est_par_err;
est_data.pss = cat(1,optimal_optsol.pss);
est_data.exitflag = exitflags;
 