% separate out integration from data processing function in combine results
function optsol = reestimate_optsol(optsol,exp_data,odep,p_id,odeopts)
% re-do perturbations in exp_data to check for consistency of solutions
% initial value for these perturbations is taken from the wt model estimate
% create parameter and options structures for these perturbations
nval = size(optsol,2); 

% combine all estimated data
est_xss = cat(1,optsol.xss);
% est_fss = cat(2,est_data.fss);
est_par = cat(2,optsol.xpar);

wt_est_xss = est_xss(:,end); % wt est data

opts = struct([]);
newodep = repmat(odep,nval,1);
newodep(:,p_id) = est_par';
for ival = 1:nval    
    opts(ival).x0 = wt_est_xss(3*(ival-1)+1:3*ival);
    opts(ival).tspan = odeopts.tspan;
    opts(ival).solver_opts = odeopts.solver_opts;
    opts(ival).odep = newodep(ival,:);    
end
npert = size(exp_data,2);
[pt_val(1:npert).exp_pid] = exp_data(1:npert).exp_pid;
[pt_val(1:npert).exp_pval] = exp_data(1:npert).exp_pval;
wt_pt_val = pt_val(end);
other_pt_val = pt_val(1:end-1);

if nval>0
    % determine ss of optimal estimated concentrations
    for ival = 1:nval
        opt_pss = ones(1,npert);        
        opt_pss(optsol(ival).xss(2,:)>optsol(ival).xss(1,:)) = 2;
        optsol(ival).pss_opt = opt_pss;
    end
    if nval>1
        parfor ival = 1:nval
            pss = ones(1,npert);
            % calculate fss for estimated xss
            odep = newodep(ival,:); 
%             wt_est_fss = kotte_flux_noCAS(wt_est_xss(3*(ival-1)+1:3*ival),odep);

            new_odeopts = opts(ival);
            % do wt perturbation first to get ss (if it already is not @ ss)

            sol_wt = getperturbations(wt_pt_val,@perturb_nonoise,new_odeopts);
            new_odeopts.x0 = sol_wt.xss;

            % recalculate perturbations for new parameters from new wt ss
            sol = getperturbations(other_pt_val,@perturb_nonoise,new_odeopts);
            close all

            % collect all solutions
            xss = [cat(2,sol.xss) sol_wt.xss];
            fss = [cat(2,sol.fss) sol_wt.fss];
            optsol(ival).xss_calc = xss;
            optsol(ival).fss_calc = fss;
            pss(xss(2,:)>xss(1,:)) = 2;
            optsol(ival).pss_calc = pss;
                        

            
        end
    else
        for ival = 1:nval
            pss = ones(1,npert);
            % calculate fss for estimated xss
            odep = newodep(ival,:); 
            wt_est_fss = kotte_flux_noCAS(wt_est_xss(3*(ival-1)+1:3*ival),odep);

            new_odeopts = opts(ival);
            % do wt perturbation first to get ss (if it already is not @ ss)

            sol_wt = getperturbations(wt_pt_val,@perturb_nonoise,new_odeopts);
            new_odeopts.x0 = sol_wt.xss;

            % recalculate perturbations for new parameters from new wt ss
            sol = getperturbations(pt_val(1:end-1),@perturb_nonoise,new_odeopts);
            close all

            % collect all solutions
            xss = [cat(2,sol.xss) sol_wt.xss];
            fss = [cat(2,sol.fss) sol_wt.fss];
            optsol(ival).xss_calc = xss;
            optsol(ival).fss_calc = fss;
            pss(xss(2,:)>xss(1,:)) = 2;
            optsol(ival).pss = pss;
        end
    end
else    
    fprintf('No optimal solution found\n');
    return
end
        


