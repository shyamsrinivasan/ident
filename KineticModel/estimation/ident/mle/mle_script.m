% get maximum likelihood estimate of all parameters for given experimental
% data (no noise)
%% load noise free data
if ~exist('C:\Users\shyam\Documents\Courses\CHE1125Project\IntegratedModels\KineticModel')
    status = 2;
    fprintf('\nLinux System\n');
else 
    status = 1;
    fprintf('\nWindows System\n');
end

if status == 1    
    load('C:/Users/shyam/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noiseless_model/pdata_nn_sep1');
elseif status == 2    
    load('~/Documents/Courses/CHE1125Project/IntegratedModels/KineticModel/estimation/noiseless_model/pdata_nn_sep1');    
end

%% collect only needed perturbations for analysis
avail_pert = size(no_noise_sol,2);
use_pert = 1;
npert = length(use_pert);
[exp_select_sol,no_noise_select_sol] = parseperturbations(no_noise_sol,use_pert);

% use wt as initial value for all perturbations
xinit = repmat(xss,npert,1);
yinit = repmat(fss,npert,1);

freq = [1:50:1500 1501:1500:3001];
pd = makedist('Normal','mu',0,'sigma',.001);
ynoise = random(pd,4,length(freq));

mle_opts = struct('nc',3,'nf',6,'npert',npert,...                  
                  'casmodelfun',@kotteCASident,...
                  'odep',odep_bkp,...
                  'tspan',opts.tspan,...
                  'freq',freq,...                  
                  'xinit',xinit,...
                  'yinit',yinit,...
                  'xexp',exp_select_sol.xdyn(:,freq),...
                  'yexp',exp_select_sol.fdyn([1 3 4 5],freq),...
                  'ynoise',ynoise);

%% get MLE of all parameters
% initial value for optimization
scale = ones(9,1);
scale(3) = 1e6;
p0 = [odep_bkp(1:5)';odep_bkp(10:13)']./scale;

MLEvals = getMLE(mle_opts,p0,scale);  
% MLE_no_noise = getMLE_CAS(mle_opts,p0,scale);